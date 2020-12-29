# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 15:39:28 2020

@author: marie
"""
import vtk
import numpy
import random

   
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
    point.SetPhiResolution(35)
    point.SetThetaResolution(35)

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
    

    points = vtk.vtkPoints()
    cellIds = vtk.vtkIdList()

    code = obbTree.IntersectWithLine(pSource, pTarget, points, cellIds)
    
 
    pointData = points.GetData()
    noPoints = pointData.GetNumberOfTuples()
    noIds = cellIds.GetNumberOfIds()
    
    assert (noPoints == noIds)
    
    pointsInter = []
    cellIdsInter = []
    for idx in range(noPoints):
        pointsInter.append(pointData.GetTuple3(idx))
        cellIdsInter.append(cellIds.GetId(idx))
    
    return pointsInter, cellIdsInter


def calcVecR(vecInc, vecNor):
    vecInc = l2n(vecInc)
    vecNor = l2n(vecNor)
    
    vecRef = vecInc - 2*numpy.dot(vecInc, vecNor)*vecNor
    
    return n2l(vecRef)

def calcAngle(A,B): 
    A = l2n(A)
    B = l2n(B)
    
    prodscal = A[0] * B[0] + A[1] * B[1] + A[2] * B[2]
    NormeA = numpy.sqrt(A[0]**2 + A[1]**2 + A[2]**2)
    NormeB = numpy.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
   
    return numpy.arccos( prodscal / (NormeA * NormeB))


def Fresnel(n1, n2, vecInc, vecNor):
    vecInc = l2n(vecInc)
    vecNor = l2n(vecNor)
    
    thetai = calcAngle(vecInc,vecNor)
    
    if thetai < numpy.arcsin((n1/n2)*sin(thetai)):
      sin_thetat_2 = ((n1/n2)**2)*numpy.sin(thetai)**2
      cos_thetat = numpy.sqrt(1-sin_thetat_2)
    
      R_parallele = ((n1*numpy.cos(thetai) - n2*cos_thetat)/(n1*numpy.cos(thetai) + n2*cos_thetat))**2
      R_perpendiculaire = ((n2*numpy.cos(thetai) - n1*cos_thetat)/(n2*numpy.cos(thetai) + n1*cos_thetat))**2
    
      R = (R_parallele + R_perpendiculaire)/2
      T = 1-R #transmis
    
    else:
       R = 1
       T = 0 
    x = random.choices(['R', 'T'], weights = [R, T])
    
    return x

def cellcenter_normal(forme):

 cellCenterCalc = vtk.vtkCellCenters()
 cellCenterCalc.SetInputConnection(forme.GetOutputPort())
 cellCenterCalc.Update()

 pointsCellCenters = cellCenterCalc.GetOutput(0)

# Create a new 'vtkPolyDataNormals' and connect to the 'lamp' half-sphere
 normalsCalc = vtk.vtkPolyDataNormals()
 normalsCalc.SetInputConnection(forme.GetOutputPort())

# Disable normal calculation at cell vertices
 normalsCalc.ComputePointNormalsOff()
# Enable normal calculation at cell centers
 normalsCalc.ComputeCellNormalsOn()
# Disable splitting of sharp edges
 normalsCalc.SplittingOff()
# Disable global flipping of normal orientation
 normalsCalc.FlipNormalsOff()
# Enable automatic determination of correct normal orientation
 normalsCalc.AutoOrientNormalsOn()
# Perform calculation
 normalsCalc.Update()
 
 normalsforme = normalsCalc.GetOutput().GetCellData().GetNormals()
    
 return pointsCellCenters, normalsforme