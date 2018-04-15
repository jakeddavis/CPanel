/*******************************************************************************
 * Copyright (c) 2014 Chris Satterwhite
 * Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
 *
 * This program and the accompanying materials are made
 * available under the terms of the Eclipse Public License 2.0
 * which is available at https://www.eclipse.org/legal/epl-2.0/
 *
 * See LICENSE.md file in the project root for full license information.
 *
 * SPDX-License-Identifier: EPL-2.0
 *
 * Contributors:
 *    Chris Satterwhite - initial code and implementation
 *    Connor Sousa - Vortex particle implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#include "wake.h"
#include "wakePanel.h"
#include "wakeLine.h"
#include "edge.h"
#include "cpNode.h"
#include "geometry.h"

#define _USE_MATH_DEFINES
#include <math.h>


wake::~wake()
{
    for (wakeLines_index_type i=0; i<wakeLines.size(); i++)
    {
        delete wakeLines[i];
    }
    for (wakePanels_index_type  i=0; i<vortexSheets.size(); i++)
    {
        delete vortexSheets[i];
    }
}

void wake::addPanel(wakePanel* wPan)
{
    Eigen::Vector3d pnt;
    std::vector<cpNode*> nodes = wPan->getNodes();
    if (wpanels.size() == 0)
    {
        // Initialize wake dimensions with first panel
        pnt = nodes[0]->getPnt();
        yMax = pnt(1);
        yMin = yMax;
        x0 = pnt(0);
        xf = x0;
        z0 = pnt(2);
        zf = z0;
        normal = wPan->getNormal();
    }
    
    for (size_t i=0; i<nodes.size(); i++)
    {
        pnt = nodes[i]->getPnt();
        if (pnt(1) > yMax)
        {
            yMax = pnt(1);
        }
        else if (pnt(1) < yMin)
        {
            yMin = pnt(1);
        }
        if (pnt(0) > xf)
        {
            xf = pnt(0);
        }
        else if (pnt(0) < x0)
        {
            x0 = pnt(0);
        }
        if (pnt(2) > zf)
        {
            zf = pnt(2);
        }
        else if (pnt(2) < z0)
        {
            z0 = pnt(2);
        }
    }
    wpanels.push_back(wPan);
}

bool wake::isSameWake(wake* other)
{
    if (other == this)
    {
        return false;
    }
    
    double eps = pow(10, -2);
    if (std::abs(other->getX0() - x0) < eps && std::abs(other->getZ0() - z0) < eps && std::abs(other->getXf() - xf) < eps && std::abs(other->getZf() - zf) < eps)
    {
        return true;
    }
    
    return false;
}

void wake::mergeWake(wake *other)
{
    std::vector<wakePanel*> pans = other->getPanels();
    wakePanel* w;
    for (size_t i=0; i<pans.size(); i++)
    {
        w = pans[i];
        wpanels.push_back(w);
        w->setParentWake(this);
    }
    
    std::vector<wakeLine*> otherLines = other->getWakeLines();
    wakeLine* wLine;
    for (size_t i=0; i<otherLines.size(); i++)
    {
        wLine = new wakeLine(*otherLines[i]);
        addWakeLine(wLine);
    }
    
    if (other->getYMin() < yMin)
    {
        yMin = other->getYMin();
    }
    if (other->getYMax() > yMax)
    {
        yMax = other->getYMax();
    }
}

void wake::addTEPanel(wakePanel* p)
{
    TEpanels.push_back(p);
}

void wake::addWakeLine(wakeLine* wl)
{
    wakeLines.push_back(wl);
    std::sort(wakeLines.begin(),wakeLines.end(),[](wakeLine* w1, wakeLine* w2) {return w1->getY() < w2->getY();});
}

void wake::trefftzPlane(double Vinf,double Sref)
{
    int nPnts = 164;
    if (nPnts % 2 != 0)
    {
        //Number is odd and needs to be even for simpsons rule integration.
        nPnts++;
    }
    Eigen::VectorXd w,dPhi;
    yLoc.resize(nPnts+1);
    yLoc(0) = yMin;
    yLoc(nPnts) = yMax;
    double step = (yMax-yMin)/(nPnts);
    w = Eigen::VectorXd::Zero(nPnts+1);
    dPhi = Eigen::VectorXd::Zero(nPnts+1);
    Cl = Eigen::VectorXd::Zero(nPnts+1);
    Cd = Eigen::VectorXd::Zero(nPnts+1);
    double xTrefftz = x0+2*(xf-x0)/3;

    
    Eigen::Vector3d pWake;
    for (int i=1; i<nPnts; i++)
    {
        yLoc(i) = yMin+i*step;
        pWake = pntInWake(xTrefftz, yLoc(i));
        w(i) = Vradial(pWake);
        dPhi(i) = -wakeStrength(yLoc(i));
        
        Cl(i) = 2*dPhi(i)/(Vinf*Sref);
        Cd(i) = dPhi(i)*w(i)/(Vinf*Vinf*Sref);
    }
    
    int i=0;
    CL = 0;
    CD = 0;
    while (i < Cl.rows()-2)
    {
        CL += 1.0/3*step*(Cl(i)+4*Cl(i+1)+Cl(i+2));
        CD += 1.0/3*step*(Cd(i)+4*Cd(i+1)+Cd(i+2));
        i += 2;
    }
    
}


void wake::trefftzPlaneVP(double Vinf,double Sref, std::vector<particle*>* particles, int numSimSteps){
    
    
    // Finding row of particles to use
    if(numSimSteps % 2 != 0){
        numSimSteps++;
    }
    int partRow = numSimSteps/2 + 2; // Particles aren't shed until timestep 3
    
    // Collect these particles in a vector
    std::vector<particle*> unsortedParts;
    for (size_t i=0; i<particles->size(); i++)
    {
        if ((*particles)[i]->shedTime == partRow)
        {
            unsortedParts.push_back((*particles)[i]);
        }
    }
    
    
    // Sort this vector from where they were shed on a y position basis, so can connect with a single continuous curve
    std::vector<particle*> sortedParts;
    while (unsortedParts.size() > 0)
    {
        // Start with first element in vector
        double lowestY = unsortedParts[0]->parentPanel->getCenter().y();
        size_t lowestYindex = 0;
        
        // Go through the vector and see if there is an element with a lesser 'y' value
        for (size_t i=0; i<unsortedParts.size(); i++) {
            double partY = unsortedParts[i]->parentPanel->getCenter().y();
            
            if(partY < lowestY){
                lowestY = partY;
                lowestYindex = i;
            }
        }
        
        sortedParts.push_back(unsortedParts[lowestYindex]);
        unsortedParts.erase(unsortedParts.begin() + static_cast<std::vector<particle *>::difference_type>(lowestYindex));
        
    }
    
    
    // Create a curve, S, to integrate upon
    
    // 1. Find curve length
    double Slen = 0;
    for(size_t i=0; i<sortedParts.size()-1; i++)
    {
        // Grab current and next particle positions
        Eigen::Vector3d pos1 = sortedParts[i]->pos;
        Eigen::Vector3d pos2 = sortedParts[i+1]->pos;
        
        // Find distance between segments and add to Slength
        Slen += (pos2 - pos1).norm();
    }
    
    
    // 2. Find points that lie on the curve. Will make 3*N particles so that there are always a odd number of pts (for simpsons integration), regardless of even or odd number of particles on S.
    
    int nPnts = (int)sortedParts.size()*3;
    double step = Slen/nPnts;
    std::vector<Eigen::Vector3d> Spts;
    std::vector<particle*> SptsP1;
    std::vector<particle*> SptsP2;
    
    Spts.push_back(sortedParts[0]->pos);
    SptsP1.push_back(sortedParts[0]);
    SptsP2.push_back(sortedParts[1]);
    
    Eigen::Vector3d nextP = sortedParts[1]->pos;
    size_t nextPindex = 1;
    Eigen::Vector3d pt = sortedParts[0]->pos;
    
    for (int i=1; i<nPnts; i++) {
        if(step < (nextP-pt).norm()){
            Eigen::Vector3d newPt = pt + step*(nextP - pt).normalized();
            Spts.push_back(newPt);
            SptsP1.push_back(sortedParts[nextPindex-1]);
            SptsP2.push_back(sortedParts[nextPindex]);
            pt = newPt;
        }
        else
        {
            double distCov = (nextP - pt).norm();
            pt = nextP; // shift pt to next node and restart
            nextP = sortedParts[nextPindex + 1]->pos;
            nextPindex ++;
            Eigen::Vector3d newPt = pt + (step-distCov)*(nextP - pt).normalized();
            Spts.push_back(newPt);
            SptsP1.push_back(sortedParts[nextPindex-1]);
            SptsP2.push_back(sortedParts[nextPindex]);
            pt = newPt;
        }
    }
    
    
    Eigen::VectorXd w,dPhi;
    w = Eigen::VectorXd::Zero(nPnts+1);
    dPhi = Eigen::VectorXd::Zero(nPnts+1);
    Cl = Eigen::VectorXd::Zero(nPnts+1);
    Cd = Eigen::VectorXd::Zero(nPnts+1);
    
    for (size_t i=1; i<Spts.size(); i++)
    {
        Eigen::Vector3d pWake = Spts[i];
        
        Eigen::Vector3d partV = Eigen::Vector3d::Zero();
        for (size_t j=0; j<particles->size(); j++) {
            partV += (*particles)[j]->velInfl(pWake);
        }
        
        double parentPanWeightedY = particlePntInWakeY( Spts[i] , SptsP1[i] , SptsP2[i] );
        double stFac = stretchFactor( SptsP1[i] , SptsP2[i] );
        Eigen::VectorXd::Index ii(static_cast<Eigen::VectorXd::Index>(i));
        
        w(ii) = std::abs(partV.z());

        dPhi(ii) = -wakeStrength(parentPanWeightedY) * stFac;
        
        Cl(ii) = 2*dPhi(ii)/(Vinf*Sref);
        Cd(ii) = dPhi(ii)*w(ii)/(Vinf*Vinf*Sref);
    }
    
    int i=0;
    CL = 0;
    CD = 0;
    while (i < Cl.rows()-2)
    {
        CL += 1.0/3*step*(Cl(i)+4*Cl(i+1)+Cl(i+2));
        CD += 1.0/3*step*(Cd(i)+4*Cd(i+1)+Cd(i+2));
        i += 2;
    }
}



double wake::particlePntInWakeY(Eigen::Vector3d pt , particle* P1 , particle* P2){
    
    //   |<-----------d----------->|
    //   |<--a-->|<---------b----->|
    //   *       *                 *
    //   P1      pt                P2
    
    
    // Find distance from pnt to P1
    double d = (P2->pos - P1->pos).norm();
    double a = (pt - P1->pos).norm();
    double b = (P2->pos - pt).norm();
    
    double smallNum = 1e-10;
    
    double ratio;
    if (a < smallNum) // pt is on P1
    {
        return P1->parentPanel->getCenter().y();
    }
    else if (b < smallNum) // pt is on P2
    {
        return P2->parentPanel->getCenter().y();
    }
    else // pt is btw P1 & P2
    {
        ratio = a/d;
    }
    
    
    // Weight y coordinate from panels
    double P1panY = P1->parentPanel->getCenter().y();
    double P2panY = P2->parentPanel->getCenter().y();
    
    return P1panY + ratio*(P2panY-P1panY);
    
}

double wake::stretchFactor(particle* P1 , particle* P2){
    
    double dShed = (P2->parentPanel->getCenter() - P1->parentPanel->getCenter()).norm();
    double dCurr = (P2->pos - P1->pos).norm();
    
    return dShed/dCurr;
    
}


double wake::dPhiWeighted(Eigen::Vector3d pt , particle* P1 , particle* P2){
    
    //   |<-----------d----------->|
    //   |<--a-->|<---------b----->|
    //   *       *                 *
    //   P1      pt                P2
    
    
    // Find distance from pnt to P1
    double d = (P2->pos - P1->pos).norm();
    double a = (pt - P1->pos).norm();
    double b = (P2->pos - pt).norm();

    
    // Find particle circulations
    double L_P1parentEdge = P1->parentPanel->edgesInOrder()[2]->length();
    double P1circ = P1->strength.y()/L_P1parentEdge;
    
    double L_P2parentEdge = P2->parentPanel->edgesInOrder()[2]->length();
    double P2circ = P2->strength.y()/L_P2parentEdge;
    
    
    // Find weighting factor
    double smallNum = 1e-10; // Rounding point error paranoia

    double ratio;
    if (a < smallNum)
    {     // pt is P1
        return P1circ;
    }
    else if (b < smallNum)
    {     // pt is P2
        return P2circ;
    }
    else
    {     // pt needs to be interpolated
        ratio = a/d;
    }
    
    return P1circ + ratio*(P2circ-P1circ);
    
}




Eigen::Vector3d wake::lambVectorInt(Eigen::VectorXd &yyLoc)
{
    // Sort by y position
    std::sort(TEpanels.begin(), TEpanels.end(), [](wakePanel* w1, wakePanel* w2) {return w1->getCenter()(1) < w2->getCenter()(1);});
    yyLoc.resize(static_cast<Eigen::VectorXd::Index>(TEpanels.size()+2));
    edge* TE = TEpanels[0]->getTE();
    
    if (TE->getN1()->getPnt()(1) > TE->getN2()->getPnt()(1))
    {
        TE->flipDir();
    }
    
    int i = 1;
    Eigen::Vector3d vel,circ;
    Eigen::MatrixXd sectForces = Eigen::MatrixXd::Zero(static_cast<Eigen::MatrixXd::Index>(TEpanels.size()+2),3);
    while (TE != nullptr)
    {
        yyLoc(i) = TE->getMidPoint()(1);
        vel = TE->edgeVelocity();
        circ = TE->TEgamma();
        sectForces.row(i) = vel.cross(circ);
        TE = TE->nextTE();
        i++;
    }
    
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    Eigen::Vector3d sectF1,sectF2;
    i = 0;
    
    while (i < sectForces.rows()-1)
    {
        sectF1 = sectForces.row(i);
        sectF2 = sectForces.row(i+1);
        F = F + 0.5*(yyLoc(i+1)-yyLoc(i))*(sectF1+sectF2);
        i++;
    }
    
    return F;
}

wakeLine* wake::findWakeLine(double y)
{
    double y1,y2;
    for (wakeLines_index_type i=0; i<wakeLines.size(); i++)
    {
        y1 = wakeLines[i]->getP1()(1);
        y2 = wakeLines[i]->getP2()(1);
        if (y <= y2 && y >= y1)
        {
            return wakeLines[i];
        }
    }
    return nullptr;
}

double wake::wakeStrength(double y)
{
    wakeLine* wl1 = nullptr;
    wakeLine* wl2 = nullptr;
    if (y < wakeLines[1]->getY())
    {
        wl1 = wakeLines[0];
        wl2 = wakeLines[1];
    }
    else if (y >= wakeLines.end()[-1]->getY())
    {
        wl1 = wakeLines.end()[-2]; //Second to last wakeline
        wl2 = wakeLines.end()[-1]; //Last wakeline
    }
    else
    {
        for (wakeLines_index_type i=1; i<wakeLines.size()-1; i++)
        {
            if ((wakeLines[i]->getY() <= y && wakeLines[i+1]->getY() > y))
            {
                wl1 = wakeLines[i];
                wl2 = wakeLines[i+1];
            }
        }
    }
    double interp = (y-wl1->getY())/(wl2->getY()-wl1->getY());
    double strength = wl1->getStrength()+interp*(wl2->getStrength()-wl1->getStrength());
    return strength;
}

double wake::Vradial(Eigen::Vector3d pWake)
{
    double r;
    double theta = M_PI/4;
    double dZmax = 0.3;
    double delZ;
    Eigen::Vector3d POI;
    POI(0) = pWake(0);
    if (pWake(1) >= 0)
    {
        r = yMax-pWake(1);
    }
    else
    {
        r = pWake(1)-yMin;
    }
    delZ = r*sin(theta);
    if (delZ > dZmax)
    {
        delZ = dZmax;
        theta = asin(dZmax/r);
    }
    if (pWake(1) >= 0)
    {
        POI(1) = yMax-r*cos(theta);
    }
    else
    {
        POI(1) = yMin+r*cos(theta);
    }
    
    POI(2) = pWake(2)+r*sin(theta);
    
    
    double Vr;
    int nPnts = 6;
    if (nPnts % 2 != 0)
    {
        nPnts++;
    }
    double dz = 0.5*delZ;
    double step = 2*dz/(nPnts-1);
    double phiPOI = 0;
    Eigen::VectorXd dPhiy(nPnts);
    Eigen::VectorXd dPhiz(nPnts);
    Eigen::MatrixXd dY(nPnts,1);
    Eigen::MatrixXd dZ(nPnts,1);
//    for (int i=0; i<wpanels.size(); i++)
//    {
//        phiPOI += wpanels[i]->panelPhi(POI);
//    }
    phiPOI = geom->wakePotential(POI);

    int i=0;
    while (i < nPnts)
    {
        double phiPnt1 = 0;
        double phiPnt2 = 0;
        double delta = -dz+i*step;
        Eigen::Vector3d ydir,zdir;
        ydir << 0,1,0;
        zdir << 0,0,1;
        Eigen::Vector3d pnt1 = POI+delta*ydir;
        Eigen::Vector3d pnt2 = POI+delta*zdir;
//        for (int j=0; j<wpanels.size(); j++)
//        {
//            phiPnt1 += wpanels[j]->panelPhi(pnt1);
//            phiPnt2 += wpanels[j]->panelPhi(pnt2);
//        }
        phiPnt1 = geom->wakePotential(pnt1);
        phiPnt2 = geom->wakePotential(pnt2);
        
        dPhiy(i) = phiPnt1-phiPOI;
        dPhiz(i) = phiPnt2-phiPOI;
        dY(i) = pnt1(1) - POI(1);
        dZ(i) = pnt2(2) - POI(2);
        i = i+1;
    }
    
    Eigen::MatrixXd Xb(0,3),Vb(0,3);
    Eigen::Vector3d V0 = Eigen::Vector3d::Zero();
    Eigen::Matrix<double,1,1> xx0;
    xx0.setZero();
    chtlsnd weightsY(xx0,dY,3,Xb,Vb,V0);
    double v = weightsY.getF().row(0)*dPhiy;
    chtlsnd weightsZ(xx0,dZ,3,Xb,Vb,V0);
    double w = weightsZ.getF().row(0)*dPhiz;
    Vr = sqrt(pow(v,2)+pow(w,2));
    return Vr;
}



Eigen::Vector3d wake::pntInWake(double x, double y)
{   //Connor's note: function finds the TE that projects out to the input point and then finds where the point lies on the trailing edge. Some testing (with a flat wake) showed that it doesn't change the x or y inputs, maybe it just finds the 'z' value?
    
    Eigen::Vector3d p1,p2,tvec,pnt,out,pntInWake;
    double t,scale;
    std::vector<edge*> edges;
    pntInWake.setZero();
    Eigen::Vector3d yDir;
    yDir << 0,1,0;
    for (wakePanels_index_type i=0; i<wpanels.size(); i++)
    {
        if (wpanels[i]->isTEpanel())
        {
            std::vector<edge*> eedges = wpanels[i]->getUpper()->getEdges();
            for (size_t j=0; j<eedges.size(); j++)
            {
                if (eedges[j]->isTE())
                {
                    p1 = eedges[j]->getNodes()[0]->getPnt();
                    p2 = eedges[j]->getNodes()[1]->getPnt();
                    if ((p1(1) <= y && p2(1) >= y) || (p1(1) >= y && p2(1) <= y))
                    {
                        t = (y-p1(1))/(p2(1)-p1(1));
                        tvec = (p2-p1);
                        pnt = p1+t*tvec;
                        out = -wpanels[i]->getNormal().cross(yDir);
                        if (out(0) < 0)
                        {
                            out *= -1; // Flip sign if p1 and p2 were out of order
                        }
                        scale = (x-pnt(0))/out(0);
                        pntInWake = pnt+scale*out;
                        return pntInWake;
                    }
                }
            }
        }
    }
    return pntInWake;
}


