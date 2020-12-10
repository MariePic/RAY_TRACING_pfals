# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 15:39:28 2020

@author: marie
"""
import vtk
import numpy

   
l2n = lambda l: numpy.array(l)
n2l = lambda n: list(n)

#### Tracer une ligner ####
def addLine(renderer, p1, p2, color=[0.0, 0.0, 1.0], opacity=1.0):
    line = vtk.vtkLineSource()
    line.SetPoint1(p1)
    line.SetPoint2(p2)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(line.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty()

    renderer.AddActor(actor)
    
#### Ajouter un point #### 
def addPoint(renderer, p, color=[0.0, 0.0, 0.0], radius=0.01, opacity=1):
    point = vtk.vtkSphereSource()
    point.SetCenter(p)
    point.SetRadius(radius)
    point.SetPhiResolution(10)
    point.SetThetaResolution(10)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(point.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)
    actor.GetProperty().EdgeVisibilityOn()  # show edges/wireframe
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty().SetEdgeColor(0, 0, 0)

    renderer.AddActor(actor)
    
    return point

### couleur lumière ###

def color(wavelenght):
    if wavelenght < 400:
        color = (0.54, 0.22, 0.69)
        
    if 400 <= wavelenght <500:
        color = (0, 0, 1)
        
    if 500 <= wavelenght <550:
        color = (0, 1, 0)
        
    if 550 <= wavelenght <600:
        color = (1, 1, 0)
        
    if 600 <= wavelenght <650:
        color = (1, 0.5, 0)
        
    if  wavelenght > 650:
        color = (1, 0, 0)
        
    return color
        
 
def isHit(obbTree, pSource, pTarget):
    r"""Returns True if the line intersects with the mesh in 'obbTree'"""
    code = obbTree.IntersectWithLine(pSource, pTarget, None, None)
    if code==0:
        return False
    return True

def GetIntersect(obbTree, pSource, pTarget):
    
    # Create an empty 'vtkPoints' object to store the intersection point coordinates
    points = vtk.vtkPoints()
    # Create an empty 'vtkIdList' object to store the ids of the cells that intersect
    # with the cast rays
    cellIds = vtk.vtkIdList()
    
    # Perform intersection
    code = obbTree.IntersectWithLine(pSource, pTarget, points, cellIds)
    
    # Get point-data 
    pointData = points.GetData()
    # Get number of intersection points found
    noPoints = pointData.GetNumberOfTuples()
    # Get number of intersected cell ids
    noIds = cellIds.GetNumberOfIds()
    
    assert (noPoints == noIds)
    
    # Loop through the found points and cells and store
    # them in lists
    pointsInter = []
    cellIdsInter = []
    for idx in range(noPoints):
        pointsInter.append(pointData.GetTuple3(idx))
        cellIdsInter.append(cellIds.GetId(idx))
    
    return pointsInter, cellIdsInter


def GetIntersect_plant(obbTree, pSource, pTarget):
    
    # Create an empty 'vtkPoints' object to store the intersection point coordinates
    points = vtk.vtkPoints()
    # Create an empty 'vtkIdList' object to store the ids of the cells that intersect
    # with the cast rays
    cellIds = vtk.vtkIdList()
    
    # Perform intersection
    code = obbTree.IntersectWithLine(pSource, pTarget, points, cellIds)
    
    # Get point-data 
    pointData = points.GetData()
    # Get number of intersection points found
    noPoints = pointData.GetNumberOfTuples()
    # Get number of intersected cell ids
    noIds = cellIds.GetNumberOfIds()
    
    assert (noPoints == noIds)
    
    # Loop through the found points and cells and store
    # them in lists
    ptsInter = []
    cIdsInter = []
    for idx in range(noPoints):
        ptsInter.append(pointData.GetTuple3(idx))
        cIdsInter.append(cellIds.GetId(idx))
    
    return ptsInter, cIdsInter


def calcVecR(vecInc, vecNor):
    vecInc = l2n(vecInc)
    vecNor = l2n(vecNor)
    
    vecRef = vecInc - 2*numpy.dot(vecInc, vecNor)*vecNor
    
    return n2l(vecRef)


def Fresnel(n1, n2, vecInc, vecNor):
    
    thetai = numpy.dot(vecInc, vecNor)
    
    sin_thetat = ((n1/n2)**2)*numpy.sin(thetai)**2
    cos_thetat = numpy.sqrt(1-(sin_thetat)**2)
    
    R_parallele = ((n1*numpy.cos(thetai) - n2*cos_thetat)/(n1*numpy.cos(thetai) + n2*cos_thetat))**2
    R_perpendiculaire = ((n2*numpy.cos(thetai) - n1*cos_thetat)/(n2*numpy.cos(thetai) + n1*cos_thetat))**2
    
    R = (R_parallele + R_perpendiculaire)/2
    T = 1-R
    
    return R, T
