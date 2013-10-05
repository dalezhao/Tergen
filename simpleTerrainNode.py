'''
Execute the following code snippets in Maya's Python script editor to try out
this node.

Please make sure that this script is already loaded by Maya's plugin manager.
'''

'''
# To generate some mountains, try:

import maya.cmds as cmds

cmds.file(new=True, f = True)
cmds.unloadPlugin("simpleTerrainNode")
cmds.loadPlugin("simpleTerrainNode")
originalMesh = cmds.polyPlane(sx = 200, sy = 200, w = 10, h = 10)[0]
simpleTerrainNode = cmds.createNode("simpleTerrainNode")
cmds.setAttr('%s.seed' % simpleTerrainNode, 0)
cmds.setAttr('%s.magnitude' % simpleTerrainNode, 8)
cmds.connectAttr('%s.outMesh' % originalMesh, '%s.inMesh' % simpleTerrainNode)
newMeshTransform = cmds.createNode("transform", n = "pTerrain")
newMesh = cmds.createNode("mesh", p = newMeshTransform)
cmds.xform(newMeshTransform, t = [0, 5, 0])
cmds.connectAttr('%s.outMesh' % simpleTerrainNode, '%s.inMesh' % newMesh)
'''

'''
# To generate a valley like terrain, try:

import maya.cmds as cmds

cmds.file(new=True, f = True)
cmds.unloadPlugin("simpleTerrainNode")
cmds.loadPlugin("simpleTerrainNode")
originalMesh = cmds.polyPlane(sx = 200, sy = 200, w = 10, h = 10)[0]
simpleTerrainNode = cmds.createNode("simpleTerrainNode")
cmds.setAttr('%s.seed' % simpleTerrainNode, 4)
cmds.setAttr('%s.magnitude' % simpleTerrainNode, 12)
cmds.setAttr('%s.lClamp' % simpleTerrainNode, -2)
cmds.setAttr('%s.lClampRange' % simpleTerrainNode, 0.5)
cmds.connectAttr('%s.outMesh' % originalMesh, '%s.inMesh' % simpleTerrainNode)
newMeshTransform = cmds.createNode("transform", n = "pTerrain")
newMesh = cmds.createNode("mesh", p = newMeshTransform)
cmds.xform(newMeshTransform, t = [0, 5, 0])
cmds.connectAttr('%s.outMesh' % simpleTerrainNode, '%s.inMesh' % newMesh)
'''

import sys
import math
import random
import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx

class SimpleTerrainNode(OpenMayaMPx.MPxNode):
    """
    A simple terrain generator dependency node.
    """
    
    # Node name and ID.
    kPluginNodeName = "simpleTerrainNode"
    kPluginNodeId = OpenMaya.MTypeId(0x00033339)
    
    # Node attribute data.

    # The input mesh based on which the terrain will be built.
    inputMeshAttr = None
    kInputMeshAttrName = "inMesh"
    kInputMeshAttrLongName = "inputMesh"
    
    # The depth of the fractal procedure used to build the terrain.
    depthAttr = None
    kDepthAttrName = "depth"
    kDepthAttrLongName = "depth"
    kDepthAttrDefault = 10
    
    # The magnitude of the terrain altitude.
    magnitudeAttr = None
    kMagnitudeAttrName = "mag"
    kMagnitudeAttrLongName = "magnitude"
    kMagnitudeAttrDefault = 5.0
    
    # The smoothness of the terrain.
    smoothAttr = None
    kSmoothAttrName = "smooth"
    kSmoothAttrLongName = "smoothness"
    kSmoothAttrDefault = 1.0
    
    # The upper bound of the altitude, above which the terrain will be clamped.
    highClampAttr = None
    kHighClampAttrName = "hClamp"
    kHighClampAttrLongName = "highClamp"
    kHighClampAttrDefault = 10.0
    
    # A range above the high clamp threshold into which the terrain above
    # the high clamp attribute will be mapped.
    highClampRangeAttr = None
    kHighClampRangeAttrName = "hClampRange"
    kHighClampRangeAttrLongName = "highClampRange"
    kHighClampRangeAttrDefault = 1.0
    
    # The lower bound of the altitude, below which the terrain will be clamped.
    lowClampAttr = None
    kLowClampAttrName = "lClamp"
    kLowClampAttrLongName = "lowClamp"
    kLowClampAttrDefault = -10.0
    
    # A range below the low clamp threshold into which the terrain below
    # the low clamp attribute will be mapped.
    lowClampRangeAttr = None
    kLowClampRangeAttrName = "lClampRange"
    kLowClampRangeAttrLongName = "lowClampRange"
    kLowClampRangeAttrDefault = 1.0
    
    # The seed value used to initialize the random generator used in the
    # terrain generation procedure.
    seedAttr = None
    kSeedAttrName = "seed"
    kSeedAttrLongName = "seed"
    kSeedAttrDefault = 0
    
    # The terrain mesh generated.
    outputMeshAttr = None
    kOutputMeshAttrName = "outMesh"
    kOutputMeshAttrLongName = "outputMesh"
    
    
    def __init__(self):
        
        OpenMayaMPx.MPxNode.__init__(self)
    
    
    @classmethod
    def nodeInitializer(cls):
        """
        Initialize this dependency node when loaded.
        """
        
        print(SimpleTerrainNode.kPluginNodeName + ": Initializing...")
        
        # Define attributes.
        
        fnMeshAttr = OpenMaya.MFnTypedAttribute()
        fnNumericAttr = OpenMaya.MFnNumericAttribute()
        
        # Input attributes.
        
        cls.inputMeshAttr = fnMeshAttr.create(
            cls.kInputMeshAttrLongName,
            cls.kInputMeshAttrName,
            OpenMaya.MFnData.kMesh
        )
        
        cls.depthAttr = fnNumericAttr.create(
            cls.kDepthAttrLongName,
            cls.kDepthAttrName,
            OpenMaya.MFnNumericData.kInt,
            cls.kDepthAttrDefault
        )
        
        cls.magnitudeAttr = fnNumericAttr.create(
            cls.kMagnitudeAttrLongName,
            cls.kMagnitudeAttrName,
            OpenMaya.MFnNumericData.kFloat,
            cls.kMagnitudeAttrDefault
        )
        fnNumericAttr.setKeyable(True)
        
        cls.smoothAttr = fnNumericAttr.create(
            cls.kSmoothAttrLongName,
            cls.kSmoothAttrName,
            OpenMaya.MFnNumericData.kFloat,
            cls.kSmoothAttrDefault
        )
        fnNumericAttr.setKeyable(True)
        
        cls.highClampAttr = fnNumericAttr.create(
            cls.kHighClampAttrLongName,
            cls.kHighClampAttrName,
            OpenMaya.MFnNumericData.kFloat,
            cls.kHighClampAttrDefault
        )
        fnNumericAttr.setKeyable(True)
        
        cls.highClampRangeAttr = fnNumericAttr.create(
            cls.kHighClampRangeAttrLongName,
            cls.kHighClampRangeAttrName,
            OpenMaya.MFnNumericData.kFloat,
            cls.kHighClampRangeAttrDefault
        )
        fnNumericAttr.setKeyable(True)
        
        cls.lowClampAttr = fnNumericAttr.create(
            cls.kLowClampAttrLongName,
            cls.kLowClampAttrName,
            OpenMaya.MFnNumericData.kFloat,
            cls.kLowClampAttrDefault
        )
        fnNumericAttr.setKeyable(True)
        
        cls.lowClampRangeAttr = fnNumericAttr.create(
            cls.kLowClampRangeAttrLongName,
            cls.kLowClampRangeAttrName,
            OpenMaya.MFnNumericData.kFloat,
            cls.kLowClampRangeAttrDefault
        )
        fnNumericAttr.setKeyable(True)
        
        cls.seedAttr = fnNumericAttr.create(
            cls.kSeedAttrLongName,
            cls.kSeedAttrName,
            OpenMaya.MFnNumericData.kInt,
            cls.kSeedAttrDefault
        )
        
        # Output attributes.
        
        cls.outputMeshAttr = fnMeshAttr.create(
            cls.kOutputMeshAttrLongName,
            cls.kOutputMeshAttrName,
            OpenMaya.MFnData.kMesh
        )
        fnMeshAttr.setWritable(False)
        fnMeshAttr.setStorable(False)
        
        # Add the attributes to the node definition.
        
        cls.addAttribute(cls.inputMeshAttr)
        cls.addAttribute(cls.depthAttr)
        cls.addAttribute(cls.magnitudeAttr)
        cls.addAttribute(cls.smoothAttr)
        cls.addAttribute(cls.highClampAttr)
        cls.addAttribute(cls.highClampRangeAttr)
        cls.addAttribute(cls.lowClampAttr)
        cls.addAttribute(cls.lowClampRangeAttr)
        cls.addAttribute(cls.seedAttr)
        cls.addAttribute(cls.outputMeshAttr)
        
        # Establishing attribute dependencies.
        
        cls.attributeAffects(cls.inputMeshAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.depthAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.magnitudeAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.smoothAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.highClampAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.highClampRangeAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.lowClampAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.lowClampRangeAttr, cls.outputMeshAttr)
        cls.attributeAffects(cls.seedAttr, cls.outputMeshAttr)
        
        print(SimpleTerrainNode.kPluginNodeName + ": Initialized.")
    
    
    @classmethod
    def nodeCreator(cls):
        
        return OpenMayaMPx.asMPxPtr(cls())
    
    
    def getTerrainCoord(self, point, bbox):
        """
        Normalize the position of a point within a bounding box volume into
        domain [0, 1] ^ 3.
        
        [params]
        
        - point: A point in the 3-D space, which should be within the bounding
                 box volume. 
        - bbox: A 3-D bounding box.
        
        [return]
        
        The normalized coordinates of the point.
        """
        
        bWidth = bbox.width()
        bHeight = bbox.height()
        bDepth = bbox.depth()
        bMinPoint = bbox.min()
        
        return [(point[0] - bMinPoint[0]) / bWidth, \
                (point[1] - bMinPoint[1]) / bHeight, \
                (point[2] - bMinPoint[2]) / bDepth]
    
    
    def genTerrainData(self, depth, magnitude, smooth, hClamp, hClampRange, lClamp, lClampRange, seed):
        """
        Compute a mesh independent terrain grid with given arguments.
        
        [params]
        
        - depth: The depth of the fractal procedure to generate the terrain.
                 The terrain grid will then have a resolution of
                 2 ^ 'depth' x 2 ^ 'depth'.
        - magnitude: The height of the terrain.
        - smooth: The smoothness of the terrain. The smaller this value, the
                  rockier the terrain.
        - hClamp: An altitude threshold above which the terrain will be clamped.
        - hClampRange: A range within which the clamped terrain above will be
                       mapped into.
        - lClamp: An altitude threshold below which the terrain will be clamped.
        - lClampRange: A range within which the clamped terrain below will be
                       mapped into.
        - seed: The seed used to initialize the random generator used in this
                procedure.
                
        [returns]
        
        The terrain grid data.
        """
        
        size0 = int(math.pow(2, depth))
        
        terrainData = [None] * (size0 + 1)
        for i in range(0, size0 + 1):
            terrainData[i] = [0] * (size0 + 1)
        
        random.seed(seed)
        
        # This is the workhorse function that recursively builds the fractal
        # terrain grid data. Each time this function will run on a grid region
        # that is half in length of the previous one, thus building details on
        # multiple levels.
        
        def genTerrainPatch(offsetu, offsetv, size):
            
            halfSize = size / 2

            # Do nothing and return if this grid region is 1 x 1 in size,
            # thus indivisible.
            if halfSize == 0:
                return
            
            # Height of the top left point.
            topLeft = terrainData[offsetv][offsetu]
            # Height of the top right point.
            topRight = terrainData[offsetv][offsetu + size]
            # Height of the bottom left point.
            bottomLeft = terrainData[offsetv + size][offsetu]
            # Height of the bottom right point.
            bottomRight = terrainData[offsetv + size][offsetu + size]
            
            # Compute the weight of the corner points used to set the middle
            # points. May change in the future to achieve different effects.
            
            def computeWeight():
                weight1 = 0.5
                weight2 = 1 - weight1
                return (weight1, weight2)
            
            # Top-middle point.
            weight1, weight2 = computeWeight()
            terrainData[offsetv][offsetu + halfSize] = topLeft * weight1 + topRight * weight2
            
            # Left-middle point.
            weight1, weight2 = computeWeight()
            terrainData[offsetv + halfSize][offsetu] = topLeft * weight1 + bottomLeft * weight2
            
            # Right-middle point.
            weight1, weight2 = computeWeight()
            terrainData[offsetv + halfSize][offsetu + size] = topRight * weight1 + bottomRight * weight2
            
            # Bottom-middle point.
            weight1, weight2 = computeWeight()
            terrainData[offsetv + size][offsetu + halfSize] = bottomLeft * weight1 + bottomRight * weight2
            
            # Middle point.
            
            # Compute the scale of the randomness of the middle point height.
            midRandomScale = magnitude * math.pow(float(halfSize) / size0, smooth)
            
            # The middle point height is computed by adding a random
            # perturbation to the average of the corner point heights.
            terrainData[offsetv + halfSize][offsetu + halfSize] = \
                (topLeft + topRight + bottomLeft + bottomRight) * 0.25 + \
                (random.random() - 0.5) * 2 * midRandomScale
            
            # Recursively calls itself until the size of the grid region is
            # 1 x 1, indivisible.
            genTerrainPatch(offsetu, offsetv, halfSize)
            genTerrainPatch(offsetu + halfSize, offsetv, halfSize)
            genTerrainPatch(offsetu, offsetv + halfSize, halfSize)
            genTerrainPatch(offsetu + halfSize, offsetv + halfSize, halfSize)
        
        # Build the terrain starting from the entire grid.
        genTerrainPatch(0, 0, size0)
        
        # Clamp the parts where heights exceed the specified clamp thresholds.
        
        for i in range(0, size0 + 1):
            
            for j in range(0, size0 + 1):
                
                height = terrainData[i][j]
                
                if height > hClamp:
                    # If the height is above the high clamp threshold, it will
                    # be normalized into the specified range above the
                    # threshold.
                    terrainData[i][j] = (1 - 1 / (abs(height - hClamp) + 1)) * hClampRange + hClamp
                    
                if height < lClamp:
                    # If the height is below the low clamp threshold, it will be
                    # normalized into the specified range below the threshold.
                    terrainData[i][j] = (1 - 1 / (abs(height - lClamp) + 1)) * lClampRange + lClamp
        
        return terrainData
    
    
    def mapTerrainData2Mesh(self, fnMesh, terrainData):
        """
        Map the data from a terrain grid onto a mesh with arbitrary vertex
        count.
        
        [params]
        
        - fnMesh: The mesh function set used to manipulate the underlying mesh.
        - terrainData: The terrain grid data, typically generated with
                       'genTerrainData' function.        
        """
        
        points = OpenMaya.MFloatPointArray()
        fnMesh.getPoints(points)
        
        # Construct a bounding box for all the vertices on the mesh.
        
        bbox = OpenMaya.MBoundingBox()
        for i in range(0, points.length()):
            point = points[i]
            bbox.expand(OpenMaya.MPoint(point))
        
        terrainSize = len(terrainData) - 1
        
        # For each point:
        
        for i in range(0, points.length()):
            
            point = points[i]
            
            # Get the normalized coordinates of the point using the bounding box
            # just constructed. Currently it is hard coded to apply the terrain
            # data along the y coordinates, thus only x and y coordinates of the
            # normalized position are used.
            tu, tv, tw = self.getTerrainCoord(point, bbox)
            
            # Get the terrain height value by linearly interpolating among the 4
            # neighbors of this point on the terrain data.
            # Also adjust for potential float point error on the boundary.
            
            ti = tu * terrainSize
            ti0 = int(ti)
            if ti0 > terrainSize:
                ti0 = terrainSize - 1
            elif ti0 < 0:
                ti0 = 0
            
            tj = tw * terrainSize
            tj0 = int(tj)
            if tj0 > terrainSize:
                tj0 = terrainSize - 1
            elif tj0 < 0:
                tj0 = 0
            
            ti1 = ti0 + 1
            tj1 = tj0 + 1
            
            # Determine the height of this point by linearly interpolating the
            # heights of its 4 neighbors on the terrain grid.
            
            # Since we assume the size of the cell on the terrain grid is 1,
            # the sum of these weights is guaranteed to be 1.
            
            c00 = (tj - tj0) * (ti - ti0)
            c01 = (tj - tj0) * (ti1 - ti)
            c10 = (tj1 - tj) * (ti - ti0)
            c11 = (tj1 - tj) * (ti1 - ti)
            
            # Points on some boundaries of the bounding box may have
            # corresponding points on the edge of the terrain grid, making some
            # of its neighbors lying out of the grid. Thus we move these
            # neighbors back on the edge of the terrain grid.
            # Note that We do it after computing the weights because otherwise
            # the weights on these edges will always be zero, making the heights
            # on these edges always zero.
            
            if ti1 > terrainSize:
                ti1 -= 1
            if tj1 > terrainSize:
                tj1 -= 1
            
            terrainValue = (c00 * terrainData[ti0][tj0] + c01 * terrainData[ti0][tj1] + \
                c10 * terrainData[ti1][tj0] + c11 * terrainData[ti1][tj1])
            
            # Offset the point by the computed height.
            #
            # Currently it is hard coded to offset the y coordinate, which will
            # be changed in the future to allow other offset directions,
            # possibly including the normal direction.
            point.y += terrainValue
            
            points.set(point, i)
        
        # Update the points on the mesh.
        fnMesh.setPoints(points)
    
    
    def compute(self, plug, dataBlock):
        
        print(SimpleTerrainNode.kPluginNodeName + ": Computing...")
        
        # Retrieves attribute values from the plugs.
                
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.depthAttr
            )
        )
        depth = dataHandle.asInt()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.magnitudeAttr
            )
        )
        magnitude = dataHandle.asFloat()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.smoothAttr
            )
        )
        smooth = dataHandle.asFloat()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.highClampAttr
            )
        )
        hClamp = dataHandle.asFloat()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.highClampRangeAttr
            )
        )
        hClampRange = dataHandle.asFloat()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.lowClampAttr
            )
        )
        lClamp = dataHandle.asFloat()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.lowClampRangeAttr
            )
        )
        lClampRange = dataHandle.asFloat()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.seedAttr
            )
        )
        seed = dataHandle.asInt()
        
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.inputValue(
                SimpleTerrainNode.inputMeshAttr
            )
        )
        inMesh = dataHandle.asMesh()
        
        # Copy the input mesh to the output plug.
        outMesh = OpenMaya.MFnMeshData().create()
        fnMesh = OpenMaya.MFnMesh(inMesh)
        fnMesh.copy(inMesh, outMesh)
        
        # Generate mesh independent terrain data.
        terrainData = self.genTerrainData(depth, magnitude, smooth, hClamp, hClampRange, lClamp, lClampRange, seed)
        
        # Map (Apply) the terrain data to the original mesh.
        self.mapTerrainData2Mesh(fnMesh, terrainData)
        
        # Write the mesh enhanced with the terrain to the output plug.
        dataHandle = OpenMaya.MDataHandle(
            dataBlock.outputValue(
                SimpleTerrainNode.outputMeshAttr
            )
        )
        dataHandle.setMObject(outMesh)
        
        # Clear the flag for this node after the computation is finished.
        dataBlock.setClean(plug)
        
        print(SimpleTerrainNode.kPluginNodeName + ": Computed.")


def initializePlugin(obj):
    
    plugin = OpenMayaMPx.MFnPlugin(
        obj,
        "Dale Zhao",
        "0.01",
        "Any"
    )
    
    try:
        
        print(SimpleTerrainNode.kPluginNodeName + ": Registering...")
        plugin.registerNode(
            SimpleTerrainNode.kPluginNodeName,
            SimpleTerrainNode.kPluginNodeId,
            SimpleTerrainNode.nodeCreator,
            SimpleTerrainNode.nodeInitializer
        )
        print(SimpleTerrainNode.kPluginNodeName + ": Registered.")
        
    except:
        
        raise Exception(SimpleTerrainNode.kPluginNodeName + ": Failed to register.")
        print(SimpleTerrainNode.kPluginNodeName + ": Failed to register.")


def uninitializePlugin(obj):
    
    plugin = OpenMayaMPx.MFnPlugin(obj)
    
    try:
        
        print(SimpleTerrainNode.kPluginNodeName + ": Deregistering...")
        plugin.deregisterNode(SimpleTerrainNode.kPluginNodeId)
        print(SimpleTerrainNode.kPluginNodeName + ": Deregistered.")
        
    except:
        
        raise Exception(SimpleTerrainNode.kPluginNodeName + ": Failed to deregister.")
        print(SimpleTerrainNode.kPluginNodeName + ": Failed to deregister.")
