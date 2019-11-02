/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
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
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#include "cpNode.h"
#include "edge.h"
#include "bodyPanel.h"


cpNode::cpNode(Eigen::Vector3d ppnt,size_t iindex) : pnt(ppnt), index(iindex), TEnode(false) {}

Eigen::Vector3d cpNode::operator-=(const cpNode &rhs)
{
    return pnt-rhs.getPnt();
}

Eigen::Vector3d operator-(cpNode lhs, const cpNode &rhs)
{
    return lhs -= rhs;
}

Eigen::Vector3d cpNode::operator+=(const cpNode &rhs)
{
    return pnt+rhs.getPnt();
}

Eigen::Vector3d operator+(cpNode lhs, const cpNode &rhs)
{
    return lhs += rhs;
}

void cpNode::addEdge(edge* e)  {edges.push_back(e);}

void cpNode::addBodyPanel(bodyPanel* p) {bodyPans.push_back(p);}

edge* cpNode::getTE(edge* exception)
{
    for (edges_index_type i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE() && edges[i] != exception)
        {
            if (edges[i]->getN2() == this)
            {
                edges[i]->flipDir();
            }
            return edges[i];
        }
    }
    return nullptr;
}

void cpNode::setTE() {TEnode = true;}

void cpNode::setIndex(size_t i) {index = i;}

edge* cpNode::getOtherTrailEdge(edge* current)
{
    for (edges_index_type i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE() && edges[i] != current)
        {
            return edges[i];
        }
    }
    return nullptr;
}


std::vector<edge*> cpNode::getTrailingEdges(){
    
    std::vector<edge*> trailingEdges;
    for (edges_index_type i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE())
        {
            trailingEdges.push_back(edges[i]);
        }
    }
    return trailingEdges;
}

Eigen::Vector3d cpNode::nodeWakeProjAngle(){
    
    std::vector<edge*> tedges = this->getTrailingEdges();
    Eigen::Vector3d avgNorm = Eigen::Vector3d::Zero();

    for (size_t j = 0; j < tedges.size(); j++) {
        std::vector<bodyPanel*> bPans = tedges[j]->getBodyPans();
        
        // In the case that a side patch is included with the top and bottom panels,
        // the displacement in the z direction will work, but not the y. CPanel seems
        // to exclude these from the node values though...
        for (bodyPanels_index_type i=0; i<bPans.size(); i++) {
            avgNorm += bPans[i]->getNormal();
        }
        
        avgNorm = avgNorm/avgNorm.size(); // Sum to average
    }
    
    return avgNorm.normalized();
}


Eigen::Vector3d cpNode::firstProjNode(double dt, double inputV){
    
    Eigen::Vector3d bisectVec = this->nodeWakeProjAngle();
    
    return this->getPnt() += c_w*dt*inputV*bisectVec;

}


Eigen::Vector3d cpNode::secProjNode(double dt, double inputV){
    
    Eigen::Vector3d bisectVec = this->nodeWakeProjAngle();

    return this->getPnt() += (c_w+1)*dt*inputV*bisectVec;

}


// lin
void cpNode::setLinCPnormal()
{
	linCPnormal.setZero();
	for (bodyPanels_index_type i = 0; i < getBodyPans().size(); i++)
	{
		linCPnormal += getBodyPans()[i]->getNormal();
	}
	linCPnormal.normalize();
}


void cpNode::setLinCPoffset()
{
	double edgeLenSum, edgeLenAvg, k;
	edgeLenSum = 0;

	for (edges_index_type i = 0; i < getEdges().size(); i++)
	{
		edgeLenSum += getEdges()[i]->length();
	}

	edgeLenAvg = edgeLenSum / getEdges().size();

	k = 0.00001; // determined from offset convergence study

	linCPoffset = k * edgeLenAvg;
}

void cpNode::setLinCPoffsetCopy()
{
	double edgeLenSum, edgeLenAvg, k;
	edgeLenSum = 0;

	for (edges_index_type i = 0; i < getEdges().size(); i++)
	{
		edgeLenSum += getEdges()[i]->length();
	}

	edgeLenAvg = edgeLenSum / getEdges().size();

	k = 0.000015; // determined from offset convergence study

	linCPoffset = k * edgeLenAvg;
}


Eigen::Vector3d cpNode::calcCP()
{
	Eigen::Vector3d ctrlPnt;
	ctrlPnt = getPnt() - linCPoffset * linCPnormal;

	return ctrlPnt;
}


void cpNode::linSetMu(double linDubStrength)
{
	linDoubletStrength = linDubStrength;
}


void cpNode::linSetPotential(Eigen::Vector3d Vinf)
{
	linPrevPotential = linPotential;

	linPotential = Vinf.dot(getPnt()) - linDoubletStrength; // Katz 11.74, 13.157
}


void cpNode::supSetPotential()
{
	linPrevPotential = linPotential;

	linPotential = linDoubletStrength / 2; // PANAIR: average potential
}
