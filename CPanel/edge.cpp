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
 *    Connor Sousa - Vortex particle implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#include "edge.h"
#include "bodyPanel.h"
#include "wakePanel.h"
#include "cpNode.h"

#define _USE_MATH_DEFINES
#include <math.h>


edge::edge(cpNode* nn1,cpNode* nn2) : n1(nn1), n2(nn2), TE(false)
{
    n1->addEdge(this);
    n2->addEdge(this);
}

void edge::addBodyPan(bodyPanel* b)
{
    
    bodyPans.push_back(b);
    checkTE();
}

void edge::addWakePan(wakePanel* w)
{
    wakePans.push_back(w);
    checkTE();
}


void edge::checkTE()
{
    if (bodyPans.size() == 2)
    {
        if (wakePans.size() == 1)
        {
            TE = true;
            bodyPanel* upper;
            bodyPanel* lower;
            Eigen::Vector3d v1,v2,normal;
            double theta1,theta2;
            normal = wakePans[0]->getNormal();
            v1 = bodyPans[0]->getCenter()-wakePans[0]->getCenter();
            v2 = bodyPans[1]->getCenter()-wakePans[0]->getCenter();
            theta1 = acos(v1.dot(normal)/(v1.norm()*normal.norm()));
            theta2 = acos(v2.dot(normal)/(v2.norm()*normal.norm()));
            
            /* By comparing the angle instead of just z location, distinction of upper and lower will be consistent even in the case of a wake shed from a vertical tail.
             *
             * \ upper
             *  \
             *   .p1   n
             *    \   /|\
             *     \___|___wake
             *     /
             *    /
             *   .p2
             *  /
             * / lower
             *
             */
            
            if (theta1 < theta2)
            {
                upper = bodyPans[0];
                lower = bodyPans[1];
            }
            else
            {
                upper = bodyPans[1];
                lower = bodyPans[0];
            }
            wakePans[0]->setParentPanels(upper,lower);
            n1->setTE();
            n2->setTE();
        }
        else
        {
            // Check for sharp edge without wake shed (i.e. vertical tail).  Used to start streamline tracing
            double angle = acos(bodyPans[0]->getNormal().dot(bodyPans[1]->getNormal()));
            if (angle > 4.7*M_PI/6 && bodyPans[0]->getID() == bodyPans[1]->getID())
//            if (angle > 5*M_PI/6 && bodyPans[0]->getID() == bodyPans[1]->getID())
            {
                TE = true;
                n1->setTE();
                n2->setTE();
                bodyPans[0]->setTEpanel(this);
                bodyPans[1]->setTEpanel(this);
            }
        }
    }
}

bool edge::isTE()  {return TE;}

bool edge::sameEdge(cpNode* node1, cpNode* node2)
{
    if ((node1 == n1 && node2 == n2) || (node1 == n2 && node2 == n1))
    {
        return true;
    }
    return false;
}

bodyPanel* edge::getOtherBodyPan(bodyPanel* currentPan)
{
    
    for (size_t i=0; i<2; i++)
    {
        if (bodyPans[i] != currentPan)
        {
            return bodyPans[i];
        }
    }
    return nullptr;
}

wakePanel* edge::getOtherWakePan(wakePanel* currentPan)
{

    for (wakePanels_index_type i=0; i<wakePans.size(); i++)
    {
        if (wakePans[i] != currentPan)
        {
            return wakePans[i];
        }
    }
    return nullptr;
}

cpNode* edge::getOtherNode(cpNode* current)
{
    if (current == n1)
    {
        return n2;
    }
    else if (current == n2)
    {
        return n1;
    }
    else
    {
        return nullptr;
    }
}

double edge::length() {return (n2->getPnt()-n1->getPnt()).norm();}

std::vector<cpNode*> edge::getNodes()
{
    std::vector<cpNode*> ns;
    ns.push_back(n1);
    ns.push_back(n2);
    return ns;
}

Eigen::Vector3d edge::getVector()
{
    return (n2->getPnt()-n1->getPnt());
}

Eigen::Vector3d edge::getMidPoint()
{
    return (n1->getPnt()+0.5*getVector());
}

Eigen::Vector3d edge::getNormal()
{
    return 0.5*(bodyPans[0]->getNormal()+bodyPans[1]->getNormal());
}

void edge::setNeighbors()
{
    if (bodyPans.size() == 2)
    {
        bodyPans[0]->addNeighbor(bodyPans[1]);
        bodyPans[1]->addNeighbor(bodyPans[0]);
    }
}

void edge::flipDir()
{
    // Flips direction of edge
    cpNode* temp;
    temp = n1;
    n1 = n2;
    n2 = temp;
}

edge* edge::nextTE()
{
    return n2->getTE(this);
}

double edge::distToEdge(const Eigen::Vector3d &pnt)
{
    Eigen::Vector3d edgeVec = n2->getPnt()-n1->getPnt();
    Eigen::Vector3d pntVec = pnt-n1->getPnt();
    double dist = (pntVec.cross(edgeVec)).norm()/edgeVec.norm();
    return dist;
}

Eigen::Vector3d edge::edgeVelocity()
{
    Eigen::Vector3d v1,v2,pnt;

    v1 = bodyPans[0]->getGlobalV();
    v2 = bodyPans[1]->getGlobalV();
    
    return 0.5*(v1+v2);
}

Eigen::Vector3d edge::TEgamma()
{
    Eigen::Vector3d gamma = Eigen::Vector3d::Zero();
    if (TE)
    {
        gamma = wakePans[0]->getMu()*getVector().normalized();
    }
    return gamma;
}

bool edge::containsNode(cpNode* node)
{
    if(node == this->getN1() || node == this->getN2())
    {
        return true;
    }
    else
    {
        return false;
    }
}



