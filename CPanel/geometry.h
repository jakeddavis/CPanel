//
//  geometry.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__geometry__
#define __CPanel__geometry__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <array>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include "panelOctree.h"
#include "surface.h"
#include "liftingSurf.h"
#include "wake.h"
#include "bodyPanel.h"
#include "wakePanel.h"
#include "edge.h"
#include "inputParams.h"
#include "cpNode.h"

class geometry
{
    std::vector<surface*> surfaces;
    std::vector<wake*> wakes;
    std::vector<bodyPanel*> bPanels;
    std::vector<wakePanel*> wPanels;
    std::vector<wakePanel*> w2Panels; // Buffer wake row two
    
    std::vector<bool> isFirstPanel;
    
    panelOctree pOctree;
    std::vector<cpNode*> nodes;
    std::vector<edge*> edges;
    short nNodes;
    short nTris;
    
    Eigen::MatrixXd A; // Doublet Influence Coefficient Matrix
    Eigen::MatrixXd B; // Source Influence Coefficient Matrix
    Eigen::MatrixXd C; // 2nd Row Buffer Wake Doublet Influence Coefficient Matrix
    
    
    bool writeCoeffFlag;
    bool vortPartFlag;
    std::string infCoeffFile;
    double dt;
    double inputV;
    
    void readTri(std::string tri_file, bool normFlag);
    std::vector<edge*> panEdges(const std::vector<cpNode*> &pNodes);
    edge* findEdge(cpNode* n1,cpNode* n2);
    void createSurfaces(const Eigen::MatrixXi &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID );
    void createOctree();
    void getLiftingSurfs(std::vector<surface*>& wakes, std::vector<surface*>& liftingSurfs);
    void setNeighbors(panel* p,int targetID);
    bool isLiftingSurf(int currentID, std::vector<int> wakeIDs);
    void correctWakeConnectivity(int wakeNodeStart,int wakeTriStart,Eigen::MatrixXi &connectivity);
    double shortestEdge(const Eigen::MatrixXi &connectivity);
    liftingSurf* getParentSurf(int wakeID);
    
    void setInfCoeff();
    Eigen::Vector4i interpIndices(std::vector<bodyPanel*> interpPans);
    
    bool infCoeffFileExists();
    void readInfCoeff();
    void writeInfCoeff();
    
    void createVPWakeSurfaces(const Eigen::MatrixXi &wakeConnectivity, const Eigen::MatrixXd &wakeNorms,  const std::vector<int> &VPwakeID, std::vector<bool> isFirstPanel);
    void calcTimeStep();
    
public:
    geometry(inputParams* p)
    {
        std::stringstream temp;
        temp << p->geomFile->name << ".infCoeff";
        infCoeffFile = temp.str();
        writeCoeffFlag = p->writeCoeffFlag;
        vortPartFlag = p->vortexParticles;
        dt = p->timeStep;
        inputV = p->velocities(0);
        
        readTri(p->geomFile->file, p->normFlag);
        
    }
    
    virtual ~geometry();
    
    geometry(const geometry& copy);
    
    geometry& operator=(const geometry &rhs);
    
    double pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf);
    double wakePotential(const Eigen::Vector3d &pnt);
    Eigen::Vector3d pntVelocity(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf, double PG);
    
    void clusterCheck();
    
    short getNumberOfNodes() {return nNodes;}
    short getNumberOfTris() {return nTris;}
    std::vector<cpNode*> getNodes() {return nodes;}
    Eigen::MatrixXd getNodePnts();
    
    std::vector<surface*> getSurfaces();
    panelOctree* getOctree() {return &pOctree;}
    std::vector<panel*> getPanels();
    std::vector<bodyPanel*>* getBodyPanels() {return &bPanels;}
    std::vector<wakePanel*>* getWakePanels() {return &wPanels;}
    std::vector<wakePanel*>* getBufferWake2Panels() {return &w2Panels;}
    std::vector<wake*> getWakes();
    Eigen::MatrixXd* getA() {return &A;}
    Eigen::MatrixXd* getB() {return &B;}
    Eigen::MatrixXd* getC() {return &C;}
    double getDt() {return dt;}
    
    void moveGeom( std::vector<double> bodyKin );
    
};


#endif /* defined(__CPanel__geometry__) */
