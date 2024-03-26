#include <vtkActor.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageData.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkVersion.h>
#include <vtkVoxelModeller.h>
#include <vtkSampleFunction.h>
#include <vtkFloatArray.h>
#include <vtkMarchingCubes.h>
#include <vtkFlyingEdges3D.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>

#include <vtkInteractorStyleTrackballActor.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkDenseArray.h>
#include <vtkCubeSource.h>
#include <vtkTriangle.h>
#include <vtkSTLWriter.h>

#include <cmath>
#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <array>
#include <eigen3/Eigen/Dense>
#include <chrono>

#include "omp.h"

using V3  = Eigen::Vector3d;
using V3f = Eigen::Vector3d;
using I3  = Eigen::Vector3i;

using nlohmann::json;

template <typename T>
T clamp(const T& n, const T& lower, const T& upper) {
  return std::max(lower, std::min(n, upper));
}

struct timer {
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
  void restart() { start_ = std::chrono::high_resolution_clock::now(); }
  timer() { restart();}
  double elapsed() {
    auto elapsed = std::chrono::high_resolution_clock::now() - start_;
    long long microseconds =
        std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    return microseconds*1.0e-6;
  }
};

struct Capsule {
  V3f start{};
  V3f end{};
  float radius{0.75};
};

void generate_images(std::vector<Capsule> &objects, const V3 &basis) {
  auto num_objects = objects.size();
  for (size_t c=0; c<num_objects; c++) {
    for (int i=-1; i<=1; i++)
      for (int j=-1; j<=1; j++)
        for (int k=-1; k<=1; k++)
          if (!(i==0 && j==0 && k==0)) {
            V3f offset{basis[0]*i, basis[1]*j, basis[2]*k};
            V3f start_pos{
                  objects[c].start.x() + offset[0],
                  objects[c].start.y() + offset[1],
                  objects[c].start.z() + offset[2]
            };
            V3f end_pos{
                  objects[c].end.x() + offset[0],
                  objects[c].end.y() + offset[1],
                  objects[c].end.z() + offset[2]
            };
            objects.push_back({start_pos, end_pos, objects[c].radius});
          }
  }
}

//Globals
std::vector<Capsule> capsules;

int main(int argc, char** argv) {
  std::cout << "FIBERGEN" << std::endl;
  std::cout << "Argc == " << argc << std::endl;
  if (argc <= 1) {
    std::cout << "Input file is required! Stopping." << std::endl;
    return 1;
  }

  std::string input_filename{argv[1]};
  std::cout << "Input file is " << input_filename << std::endl;
  std::ifstream input_file{input_filename};
  if (!input_file.is_open()) {
    std::cout << "Can't open the file. Stopping." << std::endl;
    return 1;
  }

  json input_json = json::parse(input_file);
  auto set_dvar_from_json_if_exists =
        [&input_json](const std::string field_name, double &field) {
      auto field_name_iter = input_json.find(field_name);
      if (field_name_iter != input_json.end()) field = *field_name_iter;
    };

  auto set_bvar_from_json_if_exists =
        [&input_json](const std::string field_name, bool &field) {
      auto field_name_iter = input_json.find(field_name);
      if (field_name_iter != input_json.end()) field = *field_name_iter;
    };

  auto set_ivar_from_json_if_exists =
        [&input_json](const std::string field_name, int &field) {
      auto field_name_iter = input_json.find(field_name);
      if (field_name_iter != input_json.end()) field = *field_name_iter;
    };

  bool show_viewer{false};
  set_bvar_from_json_if_exists("show_viewer", show_viewer);

  V3f cube_size = {4.0, 4.0, 4.0};
  //Reading simulation box sizes
  std::array<std::string, 3> simbox_sizes_str = {"simbox_xsize","simbox_ysize","simbox_zsize"};
  for (int i=0; i<3; i++) {
    auto simbox_size_iter = input_json.find(simbox_sizes_str[i]);
    if (simbox_size_iter == input_json.end()) {
      std::cout << "Error parsing " << simbox_sizes_str[i] << std::endl;
      return 1;
    }
    cube_size[i] = *simbox_size_iter;
  }

  std::string output_mesh_name = "output.stl";
  auto output_mesh_name_iter = input_json.find("output_mesh_name");
  if (output_mesh_name_iter != input_json.end()) {
    output_mesh_name = *output_mesh_name_iter;
  }

  auto capsule_iter = input_json.find("capsules");
  if (capsule_iter == input_json.end()) {
    std::cout << "Error parsing capsule" << std::endl;
    return 0;
  }

  auto &capsules_json = *capsule_iter;
  for (auto &capsule_json : capsules_json) {
    auto &cjs_start = capsule_json["start"];
    auto &cjs_end   = capsule_json["end"];
    float cjs_radius = float(capsule_json["radius"]);
    V3f c_start{float(cjs_start[0]), float(cjs_start[1]), float(cjs_start[2])};
    V3f c_end  {float(cjs_end[0]), float(cjs_end[1]), float(cjs_end[2])};

    std::cout << "CSX=" << c_start.x()
              << " CSY=" << c_start.y()
              << " CSZ=" << c_start.z() << std::endl;

    std::cout << "CEX=" << c_end.x()
              << " CEY=" << c_end.y()
              << " CEZ=" << c_end.z() << std::endl;

    std::cout << "CER=" << cjs_radius << std::endl << std::endl;
    capsules.push_back({c_start, c_end, cjs_radius});
  }

  bool generate_periodic{true};
  set_bvar_from_json_if_exists("periodic",  generate_periodic);
  if (generate_periodic) {
    generate_images(capsules, cube_size);
  }

  //build field
  int step_num_g = 40;
  auto step_num_iter = input_json.find("step_num");
  if (capsule_iter == input_json.end()) {
    std::cout << "Error parsing step_num" << std::endl;
    return 0;
  }
  I3 step_num  = {step_num_g, step_num_g, step_num_g};
  step_num = {int((*step_num_iter)[0]), int((*step_num_iter)[1]), int((*step_num_iter)[2])};
  V3f rstep_num = {float(step_num[0]), float(step_num[1]), float(step_num[2])};
  V3f step_size = cube_size.array() / rstep_num.array();
  std::cout << "step_size = " << step_size.transpose() << std::endl;

  int zb_steps = 4;
  set_ivar_from_json_if_exists("num_additional_steps", zb_steps);

  V3f simbox_min = V3::Zero() - zb_steps*step_size;
  V3f simbox_max = cube_size + zb_steps*step_size;
  I3 step_num_final = step_num + 2*zb_steps*I3::Ones();
  int total_steps = step_num_final[0]*step_num_final[1]*step_num_final[2];

  auto field = vtkSmartPointer<vtkImageData>::New();
  field->SetDimensions(step_num_final[0], step_num_final[1], step_num_final[2]);
  field->SetSpacing(step_size[0], step_size[1], step_size[2]);
  field->SetOrigin(simbox_min[0], simbox_min[1], simbox_min[2]);

  double mc_field_value = 1.1;
  set_dvar_from_json_if_exists("mc_field_value",  mc_field_value);

  bool clamp_field = false;
  set_bvar_from_json_if_exists("clamp_field", clamp_field);

  bool show_axes = false;
  set_bvar_from_json_if_exists("show_axes", show_axes);

  vtkSmartPointer<vtkFloatArray> fieldData = vtkSmartPointer<vtkFloatArray>::New();
  fieldData->SetNumberOfTuples(total_steps);

  timer field_build_timer;
  float OUT_OF_BOUNDS_FILL = -1111111.30;


  const int MAX_THREAD_NUM=10;
  int num_threads = 1;
  set_ivar_from_json_if_exists("num_threads",  num_threads);
  num_threads = clamp(num_threads, 1, MAX_THREAD_NUM);
  //omp_set_num_threads(num_threads);
  std::cout << "num threads is " << num_threads << std::endl;


  vtkIdType total_grid_steps = step_num_final[2]*step_num_final[2]*step_num_final[0];
  std::vector<float> field_raw_data(total_grid_steps);

  #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
  for ( int k=0; k<step_num_final[2]; k++){
    V3f pa;
    V3f ba;
    V3f pe;
    float dot_pa_ba;
    float dot_ba_ba;
    float h;
    V3f res;
    float cr;
    float field_value;
    int k_offset = k * step_num_final[1]*step_num_final[0];

    for (int j=0; j<step_num_final[1]; j++) {
      int j_offset = j * step_num_final[1];

      for (int i=0; i<step_num_final[0]; i++)  {
        //  for (int k=0; k<step_num_final[2]; k++){
        //    for (int j=0; j<step_num_final[1]; j++) {
        //      for (int i=0; i<step_num_final[0]; i++) {
        // int thread_id = omp_get_thread_num()
        int i_offset = i + k_offset + j_offset;

        if (i<zb_steps || j<zb_steps || k<zb_steps
            || (i+zb_steps > step_num_final[0])
            || (j+zb_steps > step_num_final[1])
            || (k+zb_steps > step_num_final[2])
            ) {
          field_value = OUT_OF_BOUNDS_FILL;
        } else {
          field_value  = 0.0f;

          pe[0]  =   i * step_size[0] + simbox_min[0];
          pe[1]  =   j * step_size[1] + simbox_min[1];
          pe[2]  =   k * step_size[2] + simbox_min[2];

          for (auto const &capsule : capsules) {
            pa  = pe  - capsule.start;
            ba  = capsule.end - capsule.start;
            dot_pa_ba  = pa .dot(ba );
            dot_ba_ba  = ba .dot(ba );
            h  = clamp(dot_pa_ba /dot_ba_ba , 0.0f, 1.0f);
            res  = pa  - ba *h ;
            cr  = 0.000001 + res.norm();
            field_value  += std::pow(capsule.radius/cr , 4);
          }
        }

        field_raw_data[i_offset] = field_value;
      }
    }
  }

  fieldData->SetArray(field_raw_data.data(), total_grid_steps, 1);

  omp_set_num_threads(1);

  field->GetPointData()->SetScalars(fieldData);
  std::cout << "Field built in " << field_build_timer.elapsed() << " sec." << std::endl;

  auto surface = vtkSmartPointer<vtkMarchingCubes>::New();
  surface->SetInputData(field);
  surface->ComputeNormalsOff();
  surface->SetValue(0, mc_field_value);
  surface->ReleaseDataFlagOn();
  surface->Update();

  auto colors = vtkSmartPointer<vtkNamedColors>::New();
  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());

  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->SetSize(1024, 720);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("fibergen");

  vtkNew<vtkInteractorStyleTrackballCamera> interactorStyle;
  vtkNew<vtkRenderWindowInteractor> interactor;

  interactor->SetInteractorStyle(interactorStyle);
  interactor->SetRenderWindow(renderWindow);

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(surface->GetOutputPort());
  //mapper->SetInputData(polyData);

  vtkPolyData *surface_polydata = surface->GetOutput();
  vtkPoints* surface_points = surface_polydata->GetPoints();
  double vtx[3];
  std::cout << "Number of vertices " << surface_polydata->GetNumberOfPoints() << std::endl;
  for (int i=0; i<surface_polydata->GetNumberOfPoints(); i++) {
    surface_points->GetPoint(i, vtx);
    vtx[0] = clamp(vtx[0], 0.0, cube_size[0]);
    vtx[1] = clamp(vtx[1], 0.0, cube_size[1]);
    vtx[2] = clamp(vtx[2], 0.0, cube_size[2]);
    surface_points->SetPoint(i, vtx);
    //surface_polydata->AddCellReference();
  }

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());
  actor->GetProperty()->EdgeVisibilityOn();

  renderer->AddActor(actor);

  vtkNew<vtkCubeSource> cube;
  cube->SetBounds(simbox_min[0], simbox_max[0],
                  simbox_min[1], simbox_max[1],
                  simbox_min[2], simbox_max[2]);
  cube->Update();

  vtkNew<vtkCubeSource> cube2;
  cube2->SetBounds(0.0, cube_size[0],
                   0.0, cube_size[1],
                   0.0, cube_size[2]);
  cube2->Update();

  // Mapper.
  vtkNew<vtkPolyDataMapper> cubeMapper;
  cubeMapper->SetInputData(cube->GetOutput());

  vtkNew<vtkPolyDataMapper> cubeMapper2;
  cubeMapper2->SetInputData(cube2->GetOutput());

  // Actor.
  vtkNew<vtkActor> cubeActor;
  cubeActor->SetMapper(cubeMapper);
  cubeActor->GetProperty()->EdgeVisibilityOn();
  cubeActor->GetProperty()->SetRepresentation(VTK_WIREFRAME);
  cubeActor->GetProperty()->SetColor(colors->GetColor3d("Banana").GetData());

  vtkNew<vtkActor> cubeActor2;
  cubeActor2->SetMapper(cubeMapper2);
  cubeActor2->GetProperty()->EdgeVisibilityOn();
  cubeActor2->GetProperty()->SetRepresentation(VTK_WIREFRAME);
  cubeActor2->GetProperty()->SetColor(colors->GetColor3d("Banana").GetData());

  vtkNew<vtkAxesActor> axes;

  // Assign actor to the renderer.

  if (show_axes) {
    renderer->AddActor(axes);
    renderer->AddActor(cubeActor);
  }

  renderer->AddActor(cubeActor2);

  vtkNew<vtkSTLWriter> stlWriter;
  stlWriter->SetFileName(output_mesh_name.c_str());
  stlWriter->SetInputData(surface->GetOutput(0));
  stlWriter->Write();

//  vtkNew<vtkPolyDataWriter> pdWriter;
//  pdWriter->SetFileName(output_mesh_name.c_str());
//  pdWriter->SetInputData(surface->GetOutput());
//  pdWriter->Write();

  if (show_viewer) {
    renderer->ResetCamera();
    renderWindow->Render();
    interactor->Start();
  }

  return 0;
}



