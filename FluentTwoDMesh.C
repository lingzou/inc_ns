//
//  FluentTwoDMesh.C
//
//
//  Created by Ling Zou on 9/7/13.
//
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>    // std::find
#include <sstream>      //std::ostringstream

#include "FluentTwoDMesh.h"

const static std::string TAB2  = "  ";
const static std::string TAB4  = "    ";
const static std::string TAB6  = "      ";
const static std::string TAB8  = "        ";
const static std::string TAB10 = "          ";

void FluentTriCell::addFaceAndNodes(const Face * const p_face)
{
  const long int face_id = p_face->id();
  const long int node_id1 = p_face->node_id1();
  const long int node_id2 = p_face->node_id2();

  if(std::find(_face_ids.begin(), _face_ids.end(), face_id) == _face_ids.end())
    _face_ids.push_back(face_id);
  if(std::find(_node_ids.begin(), _node_ids.end(), node_id1) == _node_ids.end())
    _node_ids.push_back(node_id1);
  if(std::find(_node_ids.begin(), _node_ids.end(), node_id2) == _node_ids.end())
    _node_ids.push_back(node_id2);
}

FluentTwoDMesh::SectionFlag
FluentTwoDMesh::extractSectionFlag(std::string line)
{
  if (line.compare(0, 2, "(0") == 0)          return COMMENTS_FLAG;
  else if (line.compare(0, 2, "(2") == 0)     return DIMENSIONS_FLAG;
  else if (line.compare(0, 3, "(10") == 0)    return NODES_FLAG;
  else if (line.compare(0, 3, "(12") == 0)    return CELLS_FLAG;
  else if (line.compare(0, 3, "(13") == 0)    return FACES_FLAG;
  // this one has to be after (10 (12 (13 etc
  else if (line.compare(0, 2, "(1") == 0)     return HEADER_FLAG;
  else if (line.compare(0, 1, "(") == 0)      return SECTIONBEGIN_FLAG;
  else if (line.compare(0, 2, "))") == 0)     return SECTIONEND_FLAG;
  else if (line.compare(0, 3, "(45") == 0)    return ZONES_FLAG;
  else                                        return UNKNOWN_FLAG;
}

void
FluentTwoDMesh::createMeshFromFile(std::string fileName, bool quiet, bool debug)
{
  if(!quiet) std::cout << "Start read Fluent mesh file: " << fileName << std::endl;

  std::vector<std::string> lines;
  std::string line;

  std::ifstream inFile;
  inFile.open(fileName.c_str());

  if(inFile.fail())
  {
    std::cerr << "ERROR: Failed to open mesh file: '" << fileName << "'\n";
    exit(1);
  }

  // read the whole file into the lines vector
  while(std::getline(inFile, line))
  {
    lines.push_back(line);
  }

  if(!quiet)
    std::cout << ">>Total line number " << lines.size() << std::endl;

  bool parsing_end = false;
  unsigned int current_line_number = 0;


  // faces are divided into different zones
  // use this face_id to give unique id to each face
  long int face_id = 0;

  while(!parsing_end)
  {
    int dummy;
    int dimension = 0;
    int N_nodes = 0,
    N_cells = 0,
    N_faces = 0;

    int type = 0;
    long int node_idx_begin = 0;
    long int node_idx_end = 0;
    int dim_node = 0;
    long int cell_idx_begin = 0;
    long int cell_idx_end = 0;
    int cell_type = 0;
    long int face_idx_begin = 0;
    long int face_idx_end = 0;
    int face_type = 0;
    int element_type = 0;

    int zone_id = 0;

    if(current_line_number > lines.size() - 1)
    {
      parsing_end = true;
      break;
    }

    std::string line = lines[current_line_number];
    if(line == "")
    {
      if(debug) std::cout << "BLANK LINE, ignored" << std::endl;
      current_line_number++;  // move to the next line
    }
    else
    {
      // std::cout << "processing line " << current_line_number + 1 << " " << line << std::endl;
      SectionFlag flg = extractSectionFlag(line);

      switch (flg)
      {
        case COMMENTS_FLAG:
          if(debug) std::cout << "COMMENTS LINE: " << line << std::endl;
          current_line_number++;  // move to the next line
          break;

        case HEADER_FLAG:
          if(debug) std::cout << "HEADER LINE: " << line << std::endl;
          current_line_number++;  // move to the next line
          break;

        case DIMENSIONS_FLAG:
          if(debug) std::cout << "DIMENSIONS LINE: " << line << std::endl;
          sscanf(line.c_str(), "(%d %d)", &dummy, &dimension);

          if(dimension != _dim)
          {
            std::cerr << "Currently deal with 2-d mesh only\n";
            exit(1);
          }

          if(debug) std::cout << "dimension is " << dimension << std::endl;
          current_line_number++;  // move to the next line
          break;

        case NODES_FLAG:
          sscanf(line.c_str(), "(10 (%d %lx %lx %d %d)", &type, &node_idx_begin, &node_idx_end, &dummy, &dim_node);
          if (type == 0)  // this is the summary line
          {
            if (node_idx_begin != 1)  // a bit sanity check
            {
              std::cerr << "ERROR: Fluent mesh, node id should start from 1. Please check your mesh file.\n";
              exit(1);
            }
            _total_Node_number = node_idx_end - node_idx_begin + 1;
            if (!quiet)
              std::cout << "NODES: (summary: " << node_idx_begin << " " << node_idx_end << ")" << std::endl;
          }
          else if (type == 1)
          {
            double x = 0, y = 0, z = 0;
            if(!quiet) std::cout << "\n>>Reading nodes data begins" << std::endl;
            bool nodes_data_finished = false;
            long int index = 0;
            while(!nodes_data_finished)
            {
              current_line_number++;
              line = lines[current_line_number];
              int n_readin = sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
              if(n_readin == dim_node)
              {
                // std::cout << "x y z " << x << " " << y << " " << z << std::endl;
                Node node(x, y, z);
                node.id() = index + 1;    // NOTE: Fluent mesh, node id starts from 1
                _NodeSet.push_back(node);
                index ++;
              }
              if(index == node_idx_end - node_idx_begin + 1)
                nodes_data_finished = true;
            }
            if(!quiet) std::cout << "<<Reading nodes data ends." << std::endl;
          }
          else
          {
            std::cerr << "I don't understand this line " << line << std::endl;
            exit(1);
          }

          current_line_number++;  // move to the next line
          break;

        case CELLS_FLAG:
          sscanf(line.c_str(), "(12 (%d %lx %lx 0 %d)", &type, &cell_idx_begin, &cell_idx_end, &cell_type);
          if (type == 0)
          {
            _total_Cell_number = cell_idx_end - cell_idx_begin + 1;
            if(!quiet)
              std::cout << "CELLS: (summary: " << cell_idx_begin << " " << cell_idx_end << ")" << std::endl;
          }
          else if (type == 1)
          {
            if(debug)
            {
              std::cout << "Cells data begins" << std::endl;
              std::cout << "We don't deal with it for 2d surface mesh." << std::endl;
            }
          }
          else
          {
            std::cerr << "I don't understand this line " << line << std::endl;
            exit(1);
          }

          current_line_number++;  // move to the next line
          break;

        case FACES_FLAG:
          sscanf(line.c_str(), "(13 (%d %lx %lx %d %d)", &zone_id, &face_idx_begin, &face_idx_end, &face_type, &element_type);
          if (zone_id == 0)
          {
            //total_Face_number = face_idx_begin - face_idx_end + 1;
            _total_Face_number = face_idx_end - face_idx_begin + 1;
            if(!quiet)
            {
              std::cout << "FACES: (summary: " << face_idx_begin << " " << face_idx_end << ")" << std::endl;
              //std::cout << "  face_idx_begin = " << face_idx_begin << std::endl;
              //std::cout << "  face_idx_end = " << face_idx_end << std::endl;
              //std::cout << "  total_Face_number = " << _total_Face_number << std::endl;
            }
          }
          else if (zone_id > 0)
          {
            //FluentFaceZone face_zone;
            //face_zone.zone_id = zone_id;
            std::vector<Face> faces;

            if(!quiet)
            {
              std::cout << "\n>>Reading faces in zone " << zone_id << std::endl;
              std::cout << "FACES: (zone: " << zone_id << " " << face_idx_begin << " " << face_idx_end << ")" << std::endl;
            }
            bool faces_data_finished = false;
            long int index = 0;
            long int node_id1, node_id2, cell_id1, cell_id2;
            while(!faces_data_finished)
            {
              current_line_number++;
              line = lines[current_line_number];
              int n_readin = sscanf(line.c_str(), "%lx %lx %lx %lx", &node_id1, &node_id2, &cell_id1, &cell_id2);
              if(n_readin == 4)
              {
                // std::cout << "n1 n2 c1 c2 " << node1 << " " << node2 << " " << cell1 << " " << cell2 << std::endl;
                Face face(node_id1, node_id2, cell_id1, cell_id2, face_id++);
                //face.id() = face_id;
                Vec3d vec = (_NodeSet[node_id2 - 1]).point() - (_NodeSet[node_id1 - 1]).point();
                face.area() = vec.norm();
                face.setFaceNormal(vec.y(), -vec.x(), 0.0);

                //FaceSet.push_back(face);
                //face_zone.Faces.push_back(face);
                faces.push_back(face);
                index ++;
                //face_id ++;
              }
              if(index == face_idx_end - face_idx_begin + 1)
                faces_data_finished = true;
            }

            //FaceZones.push_back(face_zone);
            _FaceZoneMap[zone_id] = faces;

            if(!quiet)
            {
              std::cout << "<<Reading faces in zone " << zone_id << " end" << std::endl;
            }
          }
          else
          {
            std::cerr << "I don't understand this line " << line << std::endl;
            exit(1);
          }

          current_line_number++;  // move to the next line
          break;

        case ZONES_FLAG:
          // we actually don't do anything here
          char zone_type[256];
          char zone_name[256];
          sscanf(line.c_str(), "(45 (%d %s %s)())", &zone_id, &zone_type, &zone_name);
          // std::cout << "zone_id, zone_type, zone_name " << zone_id << " " << zone_type << " " << zone_name << std::endl;
          current_line_number++;  // move to the next line
          break;

        case SECTIONBEGIN_FLAG:
        case SECTIONEND_FLAG:
          current_line_number++;  // move to the next line
          break;

        default:
          std::cerr << "I don't understand this line " << line << std::endl;
          exit(1);
      }
    }
  }

  inFile.close();

  if(!quiet) std::cout << "End of read Fluent mesh file: " << fileName << std::endl;

  /* debug */
  /*
  std::cout << "Nodes info:\n";
  for (unsigned int i = 0; i < _NodeSet.size(); i++)
    std::cout << _NodeSet[i].id() << "    " << _NodeSet[i].x() << "    " << _NodeSet[i].y() << "    " << _NodeSet[i].z() << std::endl;

  std::cout << "Faces info:\n";
  for (std::map<int, std::vector<Face> >::iterator it = _FaceZoneMap.begin(); it != _FaceZoneMap.end(); ++it)
  {
    std::cout << "Zone id: " << it->first << std::endl;
    for (unsigned int j = 0; j < it->second.size(); j++)
    {
      Face face = (it->second)[j];
      std::cout << face.id()  << "    " << face.node_id1()
                      << "    " << face.node_id2()
                      << "    " << face.cell_id1()
                      << "    " << face.cell_id2()
                      << std::endl;
    }
  }*/
  /* debug end */

  // Fluent two-d mesh data structure does not provide cell info directly.
  // Alternatively, nodes and faces are given. Adjacent cell ids of a face are given.
  // Therefore, cell data could be reconstructed while not necessarily stored.
  if(!quiet) std::cout << "\nNow processing data: " << std::endl;
  ProcessCellData();
  if(!quiet) std::cout << "Processing data end." << std::endl;

  if(!quiet) std::cout << "\nNow check if faces are properly oriented: " << std::endl;
  CheckFaceOrientation();
  if(!quiet) std::cout << "Face check end.\nEverything looks ok." << std::endl;
  if(!quiet) std::cout << "++++++++++This is the end of FluentTwoDMesh Reading++++++++++" << std::endl << std::endl;
}

void
FluentTwoDMesh::ProcessCellData()
{
  // We are going to reconstruct cell data from given faces and nodes data
  _CellSet.resize(_total_Cell_number, FluentTriCell()); //empty container

  // std::cout << "ProcessCellData start\n";
  // Loop on faces to update cell info
  for (std::map<int, std::vector<Face> >::iterator it = _FaceZoneMap.begin(); it != _FaceZoneMap.end(); ++it)
  {
    // std::cout << "Zone id: " << it->first << std::endl;
    for (unsigned int j = 0; j < (it->second).size(); j++)
    {
      Face & face = (it->second)[j];
      long int cell_id1 = face.cell_id1();
      long int cell_id2 = face.cell_id2();

      if (cell_id1 > 0)
      {
        _CellSet[cell_id1-1].addFaceAndNodes(&(it->second)[j]);
        _CellSet[cell_id1-1].addNeighborCellID(cell_id2);
      }
      if (cell_id2 > 0)
      {
        _CellSet[cell_id2-1].addFaceAndNodes(&(it->second)[j]);
        _CellSet[cell_id2-1].addNeighborCellID(cell_id1);
      }
    }
  }

  double total_volume = 0.0;
  for (long int i = 0; i < _CellSet.size(); i++)
  {
    // sanity check
    if (_CellSet[i].getNodeIDs().size() != 3)
      std::cerr << i+1 << "-th cell has nodes number: " << _CellSet[i].getNodeIDs().size() << std::endl;
    if (_CellSet[i].getFaceIDs().size() != 3)
      std::cerr << i+1 << "-th cell has faces number: " << _CellSet[i].getFaceIDs().size() << std::endl;

    // centroid
    std::vector<long int> node_ids = _CellSet[i].getNodeIDs();
    double x_centroid = 0.0, y_centroid = 0.0, z_centroid = 0.0;
    for (int k = 0; k < node_ids.size(); k++)
    {
      x_centroid += _NodeSet[node_ids[k]-1].x() / 3.0;
      y_centroid += _NodeSet[node_ids[k]-1].y() / 3.0;
      z_centroid += _NodeSet[node_ids[k]-1].z() / 3.0;
    }
    _CellSet[i].setCentroid(x_centroid, y_centroid, z_centroid);
    _CellSet[i].id() = i + 1;

    // area
    Vec3d face_vec_1 = _NodeSet[node_ids[1] - 1].point() - _NodeSet[node_ids[0] - 1].point();
    Vec3d face_vec_2 = _NodeSet[node_ids[2] - 1].point() - _NodeSet[node_ids[1] - 1].point();
    Vec3d cross_product = face_vec_1.cross(face_vec_2);
    if (cross_product.z() < 0)
      _CellSet[i].reorderNodeIDs();
    _CellSet[i].volume() = 0.5 * cross_product.norm();
    total_volume += _CellSet[i].volume();
  }
  std::cerr << "total volume = " << total_volume << std::endl;
}

void
FluentTwoDMesh::CheckFaceOrientation()
{
  for (std::map<int, std::vector<Face> >::iterator it = _FaceZoneMap.begin(); it != _FaceZoneMap.end(); ++it)
  {
    // std::cout << "Zone id: " << it->first << std::endl;
    for (unsigned int j = 0; j < (it->second).size(); j++)
    {
      Face & face = (it->second)[j];
      long int node_id1 = face.node_id1();
      long int node_id2 = face.node_id2();
      long int cell_id1 = face.cell_id1();
      long int cell_id2 = face.cell_id2();

      Vec3d face_vec = _NodeSet[node_id2 - 1].point() - _NodeSet[node_id1 - 1].point();

      if (cell_id1 > 0)
      {
        const Point ct = _CellSet[cell_id1-1].centroid();
        Vec3d node2_to_ct = ct - _NodeSet[node_id2 - 1].point();
        Vec3d cross_product = face_vec.cross(node2_to_ct);
        if (cross_product.z() < 0.0)
        {
          std::cout << "FACE ID = " << it->first << " SWAP CELLS" << std::endl;
          face.reorderCells();
        }
      }
      else
      {
        const Point ct = _CellSet[cell_id2-1].centroid();
        Vec3d node2_to_ct = ct - _NodeSet[node_id2 - 1].point();
        Vec3d cross_product = face_vec.cross(node2_to_ct);
        if (cross_product.z() > 0.0)
        {
          face.reorderCells();
          std::cout << "FACE ID = " << it->first << " SWAP CELLS" << std::endl;
        }
      }
    }
  }

  // test
  std::vector<Face> int_face = _FaceZoneMap[7];
  for (int i = 0; i < int_face.size(); i++)
  {
    Face & face = int_face[i];
    long int node_id1 = face.node_id1();
    long int node_id2 = face.node_id2();
    long int cell_id1 = face.cell_id1();
    long int cell_id2 = face.cell_id2();

    Vec3d face_vec = _NodeSet[node_id2 - 1].point() - _NodeSet[node_id1 - 1].point();
    Point ct1 = _CellSet[cell_id1-1].centroid();
    Point ct2 = _CellSet[cell_id2-1].centroid();
    Vec3d ct_to_ct = ct2 - ct1;

    Vec3d cross_product = face_vec.cross(ct_to_ct);
    if (cross_product.z() > 0.0)
      std::cout << "FACE = " << i << "WRONG." << std::endl;
  }
}

void
FluentTwoDMesh::writeMesh(FILE * ptr_File)
{
  std::ostringstream out_string_stream;

  // file head
  out_string_stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << "\n";
  out_string_stream << TAB2 << "<UnstructuredGrid>" << "\n";
  out_string_stream << TAB4 << "<Piece NumberOfPoints=\"" << _NodeSet.size() << "\" NumberOfCells=\"" << _CellSet.size() << "\">" << "\n";

  // POINTS
  out_string_stream << TAB6 << "<Points>" << "\n";    // POINTS begins

  out_string_stream << TAB8 << "<DataArray type = \"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _NodeSet.size(); i++)
    out_string_stream << TAB10 << _NodeSet[i].x() << TAB2 << _NodeSet[i].y() << TAB2 << _NodeSet[i].z() << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB6 << "</Points>" << "\n";   // POINTS ends

  // CELLS
  out_string_stream << TAB6 << "<Cells>" << "\n";     // CELLS begins

  out_string_stream << TAB8 << "<DataArray type = \"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _CellSet.size(); i++)
  {
    std::vector<long int> node_ids = _CellSet[i].getNodeIDs();
    out_string_stream << TAB10 << node_ids[0] - 1 << TAB2 << node_ids[1] - 1 << TAB2 << node_ids[2] - 1 << "\n";
  }
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB8 << "<DataArray type = \"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _CellSet.size(); i++)
    out_string_stream << TAB10 << (int)(3*i+3) << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB8 << "<DataArray type = \"UInt8\" Name=\"types\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _CellSet.size(); i++)
    out_string_stream << TAB10 << "5" << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB6 << "</Cells>" << "\n";    // CELLS ends

  fprintf(ptr_File, "%s", out_string_stream.str().c_str());
}

void
FluentTwoDMesh::finishFile(FILE * ptr_File)
{
  std::ostringstream out_string_stream;

  out_string_stream << TAB4 << "</Piece>" << "\n";
  out_string_stream << TAB2 << "</UnstructuredGrid>" << "\n";
  out_string_stream << "</VTKFile>" << "\n";

  fprintf(ptr_File, "%s", out_string_stream.str().c_str());
}

void
//FluentTwoDMesh::WriteVTUFile(FILE * ptr_File, const char* fileName)
FluentTwoDMesh::WriteVTUFile()
{
  /* FIXME */
  /* Is it potential ostringstream will be overflow? */

  // prepare file name
  std::ostringstream path, filename, fullname;
  path << "output/";
  filename << "FluentTwoDMesh.vtu";
  fullname << path.str() << filename.str();

  // write mesh info first
  FILE * ptr_File;
  ptr_File = fopen(fullname.str().c_str(), "w");








  std::ostringstream out_string_stream;

  // file head
  out_string_stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << "\n";
  out_string_stream << TAB2 << "<UnstructuredGrid>" << "\n";
  out_string_stream << TAB4 << "<Piece NumberOfPoints=\"" << _NodeSet.size() << "\" NumberOfCells=\"" << _CellSet.size() << "\">" << "\n";

  // POINTS
  out_string_stream << TAB6 << "<Points>" << "\n";    // POINTS begins

  out_string_stream << TAB8 << "<DataArray type = \"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _NodeSet.size(); i++)
    out_string_stream << TAB10 << _NodeSet[i].x() << TAB2 << _NodeSet[i].y() << TAB2 << _NodeSet[i].z() << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB6 << "</Points>" << "\n";   // POINTS ends

  // CELLS
  out_string_stream << TAB6 << "<Cells>" << "\n";     // CELLS begins

  out_string_stream << TAB8 << "<DataArray type = \"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _CellSet.size(); i++)
  {
    std::vector<long int> node_ids = _CellSet[i].getNodeIDs();
    out_string_stream << TAB10 << node_ids[0] - 1 << TAB2 << node_ids[1] - 1 << TAB2 << node_ids[2] - 1 << "\n";
  }
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB8 << "<DataArray type = \"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _CellSet.size(); i++)
    out_string_stream << TAB10 << (int)(3*i+3) << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB8 << "<DataArray type = \"UInt8\" Name=\"types\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _CellSet.size(); i++)
    out_string_stream << TAB10 << "5" << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB6 << "</Cells>" << "\n";    // CELLS ends


  // CELL DATA
//  std::ostringstream out_string_stream;
  out_string_stream << TAB6 << "<CellData>" << "\n";

  // CELL DATA (cell ID)
  out_string_stream << TAB8 << "<DataArray type=\"Float32\" Name=\"Cell_ID\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _total_Cell_number; i++)
    out_string_stream << TAB10 << _CellSet[i].id() << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";
  // CELL DATA (volume)
  out_string_stream << TAB8 << "<DataArray type=\"Float32\" Name=\"volume\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _total_Cell_number; i++)
    out_string_stream << TAB10 << _CellSet[i].volume() << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";
  // CELL DATA (temperature)
  out_string_stream << TAB8 << "<DataArray type=\"Float32\" Name=\"temperature\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _total_Cell_number; i++)
  {
    double x = _CellSet[i].centroid().x(); double y = _CellSet[i].centroid().y();
    out_string_stream << TAB10 << sin(2.0 * 3.14159265 * x) * cos(3.14159265 * y) << "\n";
  }
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB6 << "</CellData>" << "\n";


  // POINT DATA
  out_string_stream << TAB6 << "<PointData>" << "\n";
  // NODE ID
  out_string_stream << TAB8 << "<DataArray type=\"Float32\" Name=\"Node_ID\" format=\"ascii\">" << "\n";
  for(unsigned int i = 0; i < _total_Node_number; i++)
    out_string_stream << TAB10 << _NodeSet[i].id() << "\n";
  out_string_stream << TAB8 << "</DataArray>" << "\n";

  out_string_stream << TAB6 << "</PointData>" << "\n";

  // finish the file
  out_string_stream << TAB4 << "</Piece>" << "\n";
  out_string_stream << TAB2 << "</UnstructuredGrid>" << "\n";
  out_string_stream << "</VTKFile>" << "\n";

  fprintf(ptr_File, "%s", out_string_stream.str().c_str());

  fclose(ptr_File);
}
