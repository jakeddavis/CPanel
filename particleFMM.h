//
//  particleFMM.h
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 5/10/16.
//  Copyright (c) 2016 Connor Sousa. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__particleFMM__
#define __CPanel___Unstructured_Panel_Code__particleFMM__

#include <iostream>
#include "particleOctree.h"

class particleFMM
{
    // Might structure something like: build the multipole expansion, then have separate functions that either: calc velocity for a specific point or step through all of the particles and step
    
    
    particleOctree* partTree;
    short levels;
    
    void computeMultExp();

    Eigen::Vector3d findExpPos(node<particle>* thisNode);
    Eigen::Vector3d findExpPos(std::vector<particle*> parts);

    Eigen::Vector3d findExpStrength(node<particle>* thisNode);
    Eigen::Vector3d findExpStrength(std::vector<particle*> parts);

    double findExpRadius(node<particle>* thisNode);
    double findExpRadius(std::vector<particle*> parts);
    
    void printParts();
    void printMultExps();
    
    std::vector<node<particle>*> getChildren(std::vector<node<particle>*> nodeVec);


public:
    particleFMM() {};
    
    void build(particleOctree* tree);
    Eigen::Vector3d barnesHut(Eigen::Vector3d POI);
};


#endif /* defined(__CPanel___Unstructured_Panel_Code__particleFMM__) */
