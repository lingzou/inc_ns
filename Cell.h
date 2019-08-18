//
//  Cell.h
//  
//
//  Created by Ling Zou on 5/27/13.
//
//

#ifndef _Cell_h
#define _Cell_h

#include "Node.h"

class OneDCell
{
public:
  OneDCell(Node* node1, Node* node2)
  {
    node_id1 = node1->id();
    node_id2 = node2->id();
  }
  
  const unsigned int id() const {return _id;}
  const unsigned int getNodeID1() const {return node_id1;}
  const unsigned int getNodeID2() const {return node_id2;}
  
  unsigned int& id() {return _id;}
  
protected:
  unsigned int node_id1, node_id2;
  unsigned int _id;
};

#endif
