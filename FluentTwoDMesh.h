//
//  FluentTwoDMesh.h
//
//
//  Created by Ling Zou on 9/7/13.
//
//

#ifndef _FluentTwoDMesh_h
#define _FluentTwoDMesh_h

#include <string>
#include <vector>
#include <map>

#include "Node.h"
#include "Cell.h"

class Face {
public:
  Face() : _node_id1(-1), _node_id2(-1), _cell_id1(-1), _cell_id2(-1),
  _id(-1) {}
  ~Face() {}

  Face(long int node_id1, long int node_id2, long int cell_id1, long int cell_id2, long int face_id) :
    _node_id1(node_id1), _node_id2(node_id2), _cell_id1(cell_id1), _cell_id2(cell_id2),
    _id(face_id)
  {}
  const long int node_id1() const { return _node_id1; }
  const long int node_id2() const { return _node_id2; }
  const long int cell_id1() const { return _cell_id1; }
  const long int cell_id2() const { return _cell_id2; }
  const long int id() const { return _id; }
  const double area() const { return _area; }
  const Vec3d & faceNormal() const { return _face_normal; }
  inline double & area() { return _area; }
  void reorderCells()
  {
    std::swap(_cell_id1, _cell_id2);
    // flip face normal
    //setFaceNormal(-_face_normal.x(), -_face_normal.y(), 0.0);
  }
  void setFaceNormal(double x, double y, double z)
  { _face_normal.x() = x; _face_normal.y() = y; _face_normal.z() = z; }

protected:
  long int _id, _node_id1, _node_id2, _cell_id1, _cell_id2;
  double _area;
  Vec3d _face_normal;
};

class FluentTriCell {
public:
  FluentTriCell() : _id(-1), _volume(-1.0) {}
  ~FluentTriCell() {}

  inline const long int id() const { return _id; }
  inline long int & id() { return _id; }

  inline const double volume() const { return _volume; }
  inline double & volume() { return _volume; }

  inline std::vector<long int> getNodeIDs() const { return _node_ids; }
  inline std::vector<long int> getFaceIDs() const { return _face_ids; }
  inline std::vector<long int> getNbCellIDs() const { return _nb_cell_ids; }

  inline void setNodeIDs(std::vector<long int> ids) { _node_ids = ids; }
  inline void setFaceIDs(std::vector<long int> ids) { _face_ids = ids; }

  void setCentroid(double x, double y, double z) { _centroid.x() = x; _centroid.y() = y; _centroid.z() = z; }
  void setCentroid(Point point) { _centroid.x() = point.x(); _centroid.y() = point.y(); _centroid.z() = point.z(); }
  const Point & centroid() const { return _centroid; }
  void addFaceAndNodes(const Face * const p_face);
  void addNeighborCellID(long int cell_id) { _nb_cell_ids.push_back(cell_id); }
  const std::vector<long int> & getNeighborCellIDs() const { return _nb_cell_ids; }
  void reorderNodeIDs()
  {  long int tmp_id = _node_ids[2];
     _node_ids[2] = _node_ids[1];
     _node_ids[1] = tmp_id;
  }

  // Fluent two-d mesh data structure won't let us do this.
  // Alternatively, nodes and faces are given. Adjacent cell ids of a face are given.
  // Therefore, cell data could be reconstructed while not necessarily stored.
  // FluentTriCell(long int face_id1, long int face_id2, long int face_id3);

protected:
  long int _id;
  double _volume;
  Point _centroid;
  std::vector<long int> _node_ids;
  std::vector<long int> _face_ids;
  std::vector<long int> _nb_cell_ids;
};

class FluentTwoDMesh {
public:
  FluentTwoDMesh():_dim(2) {}
  ~FluentTwoDMesh() {}

  void createMeshFromFile(std::string fileName, bool quiet = true, bool debug = false);
  unsigned const int Dim() const { return _dim; }
  void ProcessCellData();
  void CheckFaceOrientation();
  // void WriteVTKFile(FILE * ptr_File, const char* fileName);
  //void WriteVTUFile(FILE * ptr_File, const char* fileName);
  void WriteVTUFile();
  void writeMesh(FILE *);
  void finishFile(FILE *);
  const unsigned int n_Nodes() const { return _total_Node_number; }
  const unsigned int n_Faces() const { return _total_Face_number; }
  const unsigned int n_Cells() const { return _total_Cell_number; }

  const std::vector<Node> & getNodeSet() const { return _NodeSet; }
  std::map<int, std::vector<Face> > & getFaceZoneMap() { return _FaceZoneMap; }
  const std::vector<FluentTriCell> & getCellSet() const { return _CellSet; }

protected:
  // File parser interfaces
  enum SectionFlag
  {
    COMMENTS_FLAG = 0,    // "(0"
    HEADER_FLAG = 1,      // "(1"
    DIMENSIONS_FLAG = 2,  // "(2"
    NODES_FLAG = 10,      // "(10"
    CELLS_FLAG = 12,      // "(12"
    FACES_FLAG = 13,      // "(13"
    SECTIONBEGIN_FLAG = 21, // "("
    SECTIONEND_FLAG = 22, // "))"
    ZONES_FLAG = 45,      // "(45"
    UNKNOWN_FLAG = 99     // unknown
  };

  SectionFlag extractSectionFlag(std::string line);

protected:
  unsigned int _dim;
  long int _total_Node_number;
  long int _total_Face_number;
  long int _total_Cell_number;

  std::vector<Node> _NodeSet;
  std::map<int, std::vector<Face> > _FaceZoneMap;
  std::vector<FluentTriCell> _CellSet;
};

#endif
