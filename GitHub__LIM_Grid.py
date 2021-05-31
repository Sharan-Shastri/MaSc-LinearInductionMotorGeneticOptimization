from GitHub__LIM_SlotPoleCalculation import *

from builtins import list

from numpy.core._multiarray_umath import ndarray

from collections import deque


class Grid(object):

    lower_slotsA: ndarray
    upper_slotsA: ndarray
    lower_slotsB: ndarray
    upper_slotsB: ndarray
    lower_slotsC: ndarray
    upper_slotsC: ndarray

    inLower_slotsA: ndarray
    outLower_slotsA: ndarray
    inUpper_slotsA: ndarray
    outUpper_slotsA: ndarray
    inLower_slotsB: ndarray
    outLower_slotsB: ndarray
    inUpper_slotsB: ndarray
    outUpper_slotsB: ndarray
    inLower_slotsC: ndarray
    outLower_slotsC: ndarray
    inUpper_slotsC: ndarray
    outUpper_slotsC: ndarray

    def __init__(self, iDesign, iSpacing, iCanvasSpacing, iMeshDensity, iN, iMeshIndexes):

        xMeshIndexes = iMeshIndexes[0]
        yMeshIndexes = iMeshIndexes[1]

        self.Spacing = iSpacing
        self.Cspacing = iCanvasSpacing
        self.meshDensity = iMeshDensity
        self.n = iN
        self.typeList = []
        self.complexTypeList = []

        self.slots = iDesign.slots
        self.poles = iDesign.poles
        self.length = iDesign.L

        if int(iDesign.slotpitch/self.Spacing) < 4:
            self.ppSlotpitch = 4
        else:
            self.ppSlotpitch = int(iDesign.slotpitch/self.Spacing)
        # Air buffer
        self.AirBuffer = iDesign.airbuffer
        if self.AirBuffer < self.Spacing:
            self.ppAirBuffer = 1
        else:
            self.ppAirBuffer = int(self.AirBuffer/self.Spacing)

        # Teeth
        self.ppTooth = int(iDesign.wt/iDesign.slotpitch*self.ppSlotpitch)
        if self.ppTooth % 2 != 0:
            self.ppTooth += 1
        # self.ppLeftEndTooth = self.ppTooth//2
        # self.ppRightEndTooth = self.ppTooth//2
        self.ppLeftEndTooth = int(iDesign.endTeeth/self.Spacing)
        self.ppRightEndTooth = self.ppLeftEndTooth

        # Slots
        self.ppSlot = self.ppSlotpitch - self.ppTooth

        # Vacuum lower
        if iDesign.vacuumBoundary < self.Spacing:
            self.ppVacuumLower = 1
        else:
            self.ppVacuumLower = int(iDesign.vacuumBoundary/self.Spacing)

        # Yoke height
        if iDesign.hy < self.Spacing:
            self.ppYokeheight = 1
        else:
            self.ppYokeheight = int(iDesign.hy/self.Spacing)

        # Slot height
        if iDesign.hs < self.Spacing:
            self.ppSlotheight = 2
        else:
            self.ppSlotheight = int(iDesign.hs/self.Spacing)
        if sum(yMeshIndexes[2]) or sum(yMeshIndexes[3]):
            if self.ppSlotheight < sum(yMeshIndexes[2])*len(self.meshDensity) + sum(yMeshIndexes[3])*len(self.meshDensity) + 2:
                self.ppSlotheight = sum(yMeshIndexes[2])*len(self.meshDensity) + sum(yMeshIndexes[3])*len(self.meshDensity) + 2
        if self.ppSlotheight % 2 != 0:
            self.ppSlotheight += 1

        # Stator height
        self.ppHeight = self.ppYokeheight + self.ppSlotheight

        # Air gap
        if iDesign.g < self.Spacing:
            self.ppAirgap = 1
        else:
            self.ppAirgap = int(iDesign.g/self.Spacing)

        # Blade rotor
        if iDesign.dr < self.Spacing:
            self.ppBladerotor = 1
        else:
            self.ppBladerotor = int(iDesign.dr/self.Spacing)

        self.Bladerotor = iDesign.dr/self.ppBladerotor
        self.BladerotorHeight = iDesign.dr/self.ppBladerotor/self.Spacing

        # Back iron
        if iDesign.bi < self.Spacing:
            self.ppBackIron = 1
        else:
            self.ppBackIron = int(iDesign.bi/self.Spacing)

        # Vacuum upper
        if iDesign.vacuumBoundary < self.Spacing:
            self.ppVacuumUpper = 1
        else:
            self.ppVacuumUpper = int(iDesign.vacuumBoundary/self.Spacing)

        self.ppLength = (self.slots - 1) * self.ppSlotpitch + self.ppSlot + self.ppLeftEndTooth + self.ppRightEndTooth + 2*self.ppAirBuffer
        self.matrix = np.array([[type('', (object,), {}) for x in np.arange(self.ppLength)] for y in np.arange(self.ppHeight + self.ppAirgap + self.ppBladerotor + self.ppBackIron + self.ppVacuumLower + self.ppVacuumUpper)])
        self.ppL = len(self.matrix[0])
        self.ppH = len(self.matrix)

        self.toothArray, self.coilArray, self.bufferArray = np.zeros((3, 1), dtype=np.int32)
        self.removeLowerCoils, self.removeUpperCoils = np.zeros((2, 1), dtype=np.int32)
        self.removeLowerCoilIdxs, self.removeUpperCoilIdxs = [], []
        self.xFirstEdgeNodes, self.xSecondEdgeNodes = [], []
        self.yFirstEdgeNodes, self.ySecondEdgeNodes = [], []

        # These lists represent whether a region boundary will have dense meshing or not. Where each region is in groups of two [#, #], one for each boundary meshing (1=use mesh density, 0=don't)
        self.xMeshSizes = np.zeros(len(xMeshIndexes), dtype=float)
        self.yMeshSizes = np.zeros(len(yMeshIndexes), dtype=float)

        # Mesh sizing
        self.fractionSize = (1 / self.meshDensity[0] + 1 / self.meshDensity[1])
        self.xListSpatialDatum = [self.AirBuffer, iDesign.endTeeth] + [iDesign.ws, iDesign.wt] * (self.slots - 1) + [iDesign.ws, iDesign.endTeeth, self.AirBuffer]
        self.xListPixelsPerRegion = [self.ppAirBuffer, self.ppLeftEndTooth] + [self.ppSlot, self.ppTooth] * (self.slots - 1) + [self.ppSlot, self.ppRightEndTooth, self.ppAirBuffer]
        self.yListSpatialDatum = [iDesign.vacuumBoundary, iDesign.hy, iDesign.hs/2, iDesign.hs/2, iDesign.g, iDesign.dr, iDesign.bi, iDesign.vacuumBoundary]
        self.yListPixelsPerRegion = [self.ppVacuumLower, self.ppYokeheight, self.ppSlotheight//2, self.ppSlotheight//2, self.ppAirgap, self.ppBladerotor, self.ppBackIron, self.ppVacuumUpper]

        for Cnt in np.arange(len(self.xMeshSizes)):
            self.xMeshSizes[Cnt] = meshBoundary([self.xListSpatialDatum[Cnt], self.xListPixelsPerRegion[Cnt]], self.Spacing, self.fractionSize, sum(xMeshIndexes[Cnt]), self.meshDensity)
            if self.xMeshSizes[Cnt] < 0:
                print('negative x mesh sizes', Cnt)

        for Cnt in np.arange(len(self.yMeshSizes)):
            self.yMeshSizes[Cnt] = meshBoundary([self.yListSpatialDatum[Cnt], self.yListPixelsPerRegion[Cnt]], self.Spacing, self.fractionSize, sum(yMeshIndexes[Cnt]), self.meshDensity)
            if self.yMeshSizes[Cnt] < 0:
                print('negative y mesh sizes')

        self.hmRegions = np.array([0, 2, 3, 4, 5], dtype=np.int32)
        self.mecRegions = np.array([1], dtype=np.int32)
        self.hmRegionsIndex = np.zeros(len(self.hmRegions) + 1, dtype=np.int32)
        self.mecRegionsIndex = np.zeros(len(self.mecRegions), dtype=np.int32)

        self.mecRegionLength = self.matrix[self.ppVacuumLower:self.ppHeight + self.ppVacuumLower, :].size

        self.vacLowerYIndexes = list(range(0, self.ppVacuumLower))
        self.vacUpperYIndexes = list(range(self.ppH-self.ppVacuumUpper, self.ppH))
        self.yokeYIndexes = list(range(self.ppVacuumLower, self.ppVacuumLower+self.ppYokeheight))
        self.upper_slotYIndexes1 = list(range(self.ppVacuumLower+self.ppYokeheight, self.ppVacuumLower+self.ppYokeheight+self.ppSlotheight//2))
        self.lower_slotYIndexes1 = list(range(self.ppVacuumLower+self.ppYokeheight+self.ppSlotheight//2, self.ppVacuumLower+self.ppYokeheight+self.ppSlotheight))
        self.airgapYIndexes = list(range(self.ppVacuumLower+self.ppHeight, self.ppVacuumLower+self.ppHeight+self.ppAirgap))
        self.bladerotorYIndexes = list(range(self.ppVacuumLower+self.ppHeight+self.ppAirgap, self.ppVacuumLower+self.ppHeight+self.ppAirgap+self.ppBladerotor))
        self.ironYIndexes = list(range(self.ppVacuumLower+self.ppHeight+self.ppAirgap+self.ppBladerotor, self.ppVacuumLower+self.ppHeight+self.ppAirgap+self.ppBladerotor+self.ppBackIron))
        self.hmYIndexes = self.vacLowerYIndexes + self.vacUpperYIndexes + self.airgapYIndexes + self.bladerotorYIndexes + self.ironYIndexes
        self.mecYIndexes = self.yokeYIndexes + self.upper_slotYIndexes1 + self.lower_slotYIndexes1

        # Thrust of the entire integration region
        self.Fx = 0.0
        self.Fy = 0.0

    def produceGrid(self, iDesign, grid, pixelSpacing, iMeshIndexes):

        xMeshIndexes = iMeshIndexes[0]
        yMeshIndexes = iMeshIndexes[1]

        #  Initialize grid with Nodes
        a, b = 0, 0
        slotCount = 0
        listOffset = self.ppAirBuffer + self.ppLeftEndTooth
        oldTooth, newTooth, oldSlot, newSlot = 0, listOffset, 0, 0
        fullToothArray, slotArray = [], []
        while a < self.ppH:
            while b < self.ppL:
                # Create indexes for slots and teeth
                if slotCount < self.slots and listOffset < b < self.ppL-listOffset:
                    oldSlot = newTooth
                    newSlot = oldSlot + self.ppSlot
                    oldTooth = newSlot
                    newTooth = oldTooth + self.ppTooth
                    if slotCount < self.slots - 1:
                        fullToothArray += list(range(oldTooth, newTooth))
                    slotCount += 1
                else:
                    pass
                b += 1
            b = 0
            a += 1

        leftEndTooth = list(range(self.ppAirBuffer, listOffset))
        rightEndTooth = list(range(fullToothArray[-1] + self.ppSlot + 1, fullToothArray[-1] + self.ppSlot + self.ppRightEndTooth + 1))
        self.toothArray = leftEndTooth + fullToothArray + rightEndTooth
        self.coilArray = [x for x in range(grid.ppAirBuffer, grid.ppL - grid.ppAirBuffer) if x not in self.toothArray]
        self.bufferArray = [x for x in list(range(grid.ppL)) if x not in self.coilArray and x not in self.toothArray]

        offset = self.ppSlot*math.ceil(iDesign.q)
        intSlots = math.ceil(iDesign.q)*iDesign.poles*3

        # Split slots into respective phases
        windingShift = 2
        temp_upper_slotArray = self.coilArray
        upper_slotArray = temp_upper_slotArray + [None] * self.ppSlot * (intSlots - self.slots)
        lower_slotArray = deque(upper_slotArray)
        lower_slotArray.rotate(-windingShift*offset)
        lower_slotArray = list(lower_slotArray)

        upper_slotArrayA, upper_slotArrayB, upper_slotArrayC = [], [], []
        lower_slotArrayA, lower_slotArrayB, lower_slotArrayC = [], [], []
        for threeSlots in np.arange(0, self.slots, 3):
            upper_slotArrayA += upper_slotArray[(threeSlots+1)*offset:(threeSlots+2)*offset]
            upper_slotArrayB += upper_slotArray[threeSlots*offset:(threeSlots+1)*offset]
            upper_slotArrayC += upper_slotArray[(threeSlots+2)*offset:(threeSlots+3)*offset]

            lower_slotArrayA += lower_slotArray[(threeSlots+1)*offset:(threeSlots+2)*offset]
            lower_slotArrayB += lower_slotArray[threeSlots*offset:(threeSlots+1)*offset]
            lower_slotArrayC += lower_slotArray[(threeSlots+2)*offset:(threeSlots+3)*offset]

        self.removeLowerCoils = [0, 1, 2, self.slots-1]
        self.removeUpperCoils = [0, self.slots-1, self.slots-2, self.slots-3]

        for idx in self.removeLowerCoils:
            coilOffset = idx*self.ppSlot
            self.removeLowerCoilIdxs += self.coilArray[coilOffset:coilOffset+self.ppSlot]

        for idx in self.removeUpperCoils:
            coilOffset = idx*self.ppSlot
            self.removeUpperCoilIdxs += self.coilArray[coilOffset:coilOffset+self.ppSlot]

        upper_slotArrayA = [x for x in upper_slotArrayA if x not in self.removeUpperCoilIdxs]
        upper_slotArrayB = [x for x in upper_slotArrayB if x not in self.removeUpperCoilIdxs]
        upper_slotArrayC = [x for x in upper_slotArrayC if x not in self.removeUpperCoilIdxs]
        lower_slotArrayA = [x for x in lower_slotArrayA if x not in self.removeLowerCoilIdxs]
        lower_slotArrayB = [x for x in lower_slotArrayB if x not in self.removeLowerCoilIdxs]
        lower_slotArrayC = [x for x in lower_slotArrayC if x not in self.removeLowerCoilIdxs]

        self.upper_slotsA = np.array(upper_slotArrayA)
        self.lower_slotsA = np.array(lower_slotArrayA)

        self.upper_slotsB = np.array(upper_slotArrayB)
        self.lower_slotsB = np.array(lower_slotArrayB)

        self.lower_slotsC = np.array(lower_slotArrayC)
        self.upper_slotsC = np.array(upper_slotArrayC)

        # Remove None from
        lowerCoilsA_NoneRemoved = list(filter(None, self.lower_slotsA))
        lowerCoilsA = np.split(np.array(lowerCoilsA_NoneRemoved), len(lowerCoilsA_NoneRemoved)//self.ppSlot, axis=0)
        upperCoilsA_NoneRemoved = list(filter(None, self.upper_slotsA))
        upperCoilsA = np.split(np.array(upperCoilsA_NoneRemoved), len(upperCoilsA_NoneRemoved)//self.ppSlot, axis=0)

        lowerCoilsB_NoneRemoved = list(filter(None, self.lower_slotsB))
        lowerCoilsB = np.split(np.array(lowerCoilsB_NoneRemoved), len(lowerCoilsB_NoneRemoved)//self.ppSlot, axis=0)
        upperCoilsB_NoneRemoved = list(filter(None, self.upper_slotsB))
        upperCoilsB = np.split(np.array(upperCoilsB_NoneRemoved), len(upperCoilsB_NoneRemoved)//self.ppSlot, axis=0)

        lowerCoilsC_NoneRemoved = list(filter(None, self.lower_slotsC))
        lowerCoilsC = np.split(np.array(lowerCoilsC_NoneRemoved), len(lowerCoilsC_NoneRemoved)//self.ppSlot, axis=0)
        upperCoilsC_NoneRemoved = list(filter(None, self.upper_slotsC))
        upperCoilsC = np.split(np.array(upperCoilsC_NoneRemoved), len(upperCoilsC_NoneRemoved)//self.ppSlot, axis=0)

        # Sort coils into direction of current (ex, in and out of page)
        self.inLower_slotsA = np.array(lowerCoilsA[1::2])
        self.outLower_slotsA = np.array(lowerCoilsA[::2])
        self.inUpper_slotsA = np.array(upperCoilsA[1::2])
        self.outUpper_slotsA = np.array(upperCoilsA[::2])
        self.inLower_slotsB = np.array(lowerCoilsB[1::2])
        self.outLower_slotsB = np.array(lowerCoilsB[::2])
        self.inUpper_slotsB = np.array(upperCoilsB[1::2])
        self.outUpper_slotsB = np.array(upperCoilsB[::2])
        self.inLower_slotsC = np.array(lowerCoilsC[::2])
        self.outLower_slotsC = np.array(lowerCoilsC[1::2])
        self.inUpper_slotsC = np.array(upperCoilsC[::2])
        self.outUpper_slotsC = np.array(upperCoilsC[1::2])

        # X Mesh Density
        Cnt = 0
        idxOffset = self.xListPixelsPerRegion[Cnt]
        idxLeft, idxRight = 0, idxOffset
        idxList = range(grid.ppL)
        # iXmeshIndexes is a list of 1s and 0s for each boundary for each region to say whether or not dense meshing is required at the boundary
        for boundary in xMeshIndexes:
            if boundary[0]:  # Left Boundary in the region
                firstIndexes = idxList[idxLeft]
                secondIndexes = idxList[idxLeft + 1]
                self.xFirstEdgeNodes.append(firstIndexes)
                self.xSecondEdgeNodes.append(secondIndexes)
            if boundary[1]:  # Right Boundary in the region
                firstIndexes = idxList[idxRight - 1]
                secondIndexes = idxList[idxRight - 2]
                self.xFirstEdgeNodes.append(firstIndexes)
                self.xSecondEdgeNodes.append(secondIndexes)
            idxLeft += idxOffset
            if Cnt < len(self.xListPixelsPerRegion) - 1:
                idxOffset = self.xListPixelsPerRegion[Cnt+1]
            idxRight += idxOffset
            Cnt += 1

        # Z Mesh Density
        Cnt = 0
        idxOffset = self.yListPixelsPerRegion[Cnt]
        idxLeft, idxRight = 0, idxOffset
        idxList = range(grid.ppH)
        for boundary in yMeshIndexes:
            if boundary[0]:  # Left Boundary in the region
                firstIndexes = idxList[idxLeft]
                secondIndexes = idxList[idxLeft + 1]
                self.yFirstEdgeNodes.append(firstIndexes)
                self.ySecondEdgeNodes.append(secondIndexes)
            if boundary[1]:  # Right Boundary in the region
                firstIndexes = idxList[idxRight - 1]
                secondIndexes = idxList[idxRight - 2]
                self.yFirstEdgeNodes.append(firstIndexes)
                self.ySecondEdgeNodes.append(secondIndexes)
            idxLeft += idxOffset
            if Cnt < len(self.yListPixelsPerRegion) - 1:
                idxOffset = self.yListPixelsPerRegion[Cnt+1]
            idxRight += idxOffset
            Cnt += 1

        xBoundaryList = [self.bufferArray[self.ppAirBuffer - 1], self.toothArray[self.ppLeftEndTooth - 1]] + self.coilArray[self.ppSlot-1::self.ppSlot] + self.toothArray[self.ppLeftEndTooth + self.ppTooth - 1:-self.ppRightEndTooth:self.ppTooth] + [self.toothArray[-1], self.bufferArray[-1]]
        yVac1Boundary = self.ppVacuumLower-1
        yYokeBoundary = yVac1Boundary + self.ppYokeheight
        yLowerCoilsBoundary = yYokeBoundary + self.ppSlotheight//2
        yUpperCoilsBoundary = yLowerCoilsBoundary + self.ppSlotheight//2
        yAirBoundary = yUpperCoilsBoundary + self.ppAirgap
        yBladeBoundary = yAirBoundary + self.ppBladerotor
        yBackIronBoundary = yBladeBoundary + self.ppBackIron
        yVac2Boundary = yBackIronBoundary + self.ppVacuumUpper
        yBoundaryList = [yVac1Boundary, yYokeBoundary, yLowerCoilsBoundary, yUpperCoilsBoundary, yAirBoundary, yBladeBoundary, yBackIronBoundary, yVac2Boundary]

        a, b = 0, 0
        c, d, e, f = 0, 0, 0, 0
        Cnt = 0
        yCnt = 0

        # Assign spatial data to the nodes
        while a < self.ppH:

            xCnt = 0
            # Keep track of the y coordinate for each node
            if a in self.yFirstEdgeNodes:
                delZ = pixelSpacing / self.meshDensity[0]
            elif a in self.ySecondEdgeNodes:
                delZ = pixelSpacing / self.meshDensity[1]
            else:
                # Vacuum lower
                if a in self.vacLowerYIndexes:
                    delZ = self.yMeshSizes[d]
                # Yoke
                elif a in self.yokeYIndexes:
                    delZ = self.yMeshSizes[d]
                # Lower coils
                elif a in self.lower_slotYIndexes1:
                    delZ = self.yMeshSizes[d]
                # Upper coils
                elif a in self.upper_slotYIndexes1:
                    delZ = self.yMeshSizes[d]
                # Air gap
                elif a in self.airgapYIndexes:
                    delZ = self.yMeshSizes[d]
                # Blade rotor
                elif a in self.bladerotorYIndexes:
                    delZ = self.yMeshSizes[d]
                # Back iron
                elif a in self.ironYIndexes:
                    delZ = self.yMeshSizes[d]
                # Vacuum upper
                elif a in self.vacUpperYIndexes:
                    delZ = self.yMeshSizes[d]
                else:
                    delZ = pixelSpacing

            while b < self.ppL:

                # Keep track of the x coordinate for each node
                if b in self.xFirstEdgeNodes:
                    delX = pixelSpacing / self.meshDensity[0]
                elif b in self.xSecondEdgeNodes:
                    delX = pixelSpacing / self.meshDensity[1]
                else:
                    # Air buffer
                    if b in self.bufferArray:
                        delX = self.xMeshSizes[c]
                    # End teeth
                    elif b in self.toothArray[:self.ppLeftEndTooth] + self.toothArray[-self.ppRightEndTooth:]:
                        delX = self.xMeshSizes[c]
                    # Coils
                    elif b in self.coilArray:
                        delX = self.xMeshSizes[c]
                    # Full teeth
                    elif b in self.toothArray and b not in self.toothArray[:self.ppLeftEndTooth] + self.toothArray[-self.ppRightEndTooth:]:
                        delX = self.xMeshSizes[c]
                    else:
                        delX = pixelSpacing

                # delX can be negative which means x can be 0
                self.matrix[a][b] = Node([b, a], [xCnt, delX], [yCnt, delZ], iDesign)

                # Keep track of the x coordinate for each node
                if b in self.xFirstEdgeNodes:
                    xCnt += pixelSpacing / self.meshDensity[0]
                elif b in self.xSecondEdgeNodes:
                    xCnt += pixelSpacing / self.meshDensity[1]
                else:
                    # Air buffer
                    if b in self.bufferArray:
                        xCnt += self.xMeshSizes[e]
                    # End teeth
                    elif b in self.toothArray[:self.ppLeftEndTooth] + self.toothArray[-self.ppRightEndTooth:]:
                        xCnt += self.xMeshSizes[e]
                    # Coils
                    elif b in self.coilArray:
                        xCnt += self.xMeshSizes[e]
                    # Full teeth
                    elif b in self.toothArray and b not in self.toothArray[:self.ppLeftEndTooth] + self.toothArray[-self.ppRightEndTooth:]:
                        xCnt += self.xMeshSizes[e]
                    else:
                        xCnt += pixelSpacing

                if b in xBoundaryList:
                    c += 1
                    e += 1

                b += 1
                Cnt += 1
            c, e = 0, 0

            # Keep track of the y coordinate for each node
            if a in self.yFirstEdgeNodes:
                yCnt += pixelSpacing / self.meshDensity[0]
            elif a in self.ySecondEdgeNodes:
                yCnt += pixelSpacing / self.meshDensity[1]
            else:
                # Vacuum lower
                if a in self.vacLowerYIndexes:
                    yCnt += self.yMeshSizes[f]
                # Yoke
                elif a in self.yokeYIndexes:
                    yCnt += self.yMeshSizes[f]
                # Lower coils
                elif a in self.lower_slotYIndexes1:
                    yCnt += self.yMeshSizes[f]
                # Upper coils
                elif a in self.upper_slotYIndexes1:
                    yCnt += self.yMeshSizes[f]
                # Air gap
                elif a in self.airgapYIndexes:
                    yCnt += self.yMeshSizes[f]
                # Blade rotor
                elif a in self.bladerotorYIndexes:
                    yCnt += self.yMeshSizes[f]
                # Back iron
                elif a in self.ironYIndexes:
                    yCnt += self.yMeshSizes[f]
                # Vacuum upper
                elif a in self.vacUpperYIndexes:
                    yCnt += self.yMeshSizes[f]
                else:
                    yCnt += pixelSpacing

            if a in yBoundaryList:
                d += 1
                f += 1

            b = 0
            a += 1

        # Assign property data to the nodes
        a, b = 0, 0
        while a < self.ppH:
            while b < self.ppL:
                if a in self.yokeYIndexes and b not in self.bufferArray:
                    self.matrix[a][b].material = 'iron'
                    self.matrix[a][b].ur = 1000.0
                    self.matrix[a][b].sigma = 4.5*10**6
                elif a in self.airgapYIndexes:
                    self.matrix[a][b].material = 'vacuum'
                    self.matrix[a][b].ur = 1.0
                    self.matrix[a][b].sigma = 3*10**(-15)
                elif a in self.vacLowerYIndexes or a in self.vacUpperYIndexes:
                    self.matrix[a][b].material = 'vacuum'
                    self.matrix[a][b].ur = 1.0
                    self.matrix[a][b].sigma = 3*10**(-15)
                elif a in self.bladerotorYIndexes:
                    self.matrix[a][b].material = 'aluminum'
                    self.matrix[a][b].ur = 1.000021
                    # https: // www.tibtech.com / conductivite.php?lang = en_US
                    self.matrix[a][b].sigma = 17.0*10**6
                elif a in self.ironYIndexes:
                    self.matrix[a][b].material = 'iron'
                    self.matrix[a][b].ur = 1000.0
                    # https: // www.thoughtco.com / table - of - electrical - resistivity - conductivity - 608499
                    self.matrix[a][b].sigma = 4.5*10**6
                else:
                    if a in self.upper_slotYIndexes1:
                        aIdx = self.upper_slotsA
                        bIdx = self.upper_slotsB
                        cIdx = self.upper_slotsC
                    elif a in self.lower_slotYIndexes1:
                        aIdx = self.lower_slotsA
                        bIdx = self.lower_slotsB
                        cIdx = self.lower_slotsC
                    else:
                        aIdx = []
                        bIdx = []
                        cIdx = []

                    if b in self.toothArray:
                        self.matrix[a][b].material = 'iron'
                        self.matrix[a][b].ur = 1000.0
                        self.matrix[a][b].sigma = 4.5*10**6
                    elif b in self.bufferArray:
                        self.matrix[a][b].material = 'vacuum'
                        self.matrix[a][b].ur = 1.0
                        self.matrix[a][b].sigma = 3*10**(-15)
                    elif b in aIdx:
                        self.matrix[a][b].material = 'copperA'
                        self.matrix[a][b].ur = 0.999991
                        self.matrix[a][b].sigma = 5.96*10**6
                    elif b in bIdx:
                        self.matrix[a][b].material = 'copperB'
                        self.matrix[a][b].ur = 0.999991
                        self.matrix[a][b].sigma = 5.96*10**6
                    elif b in cIdx:
                        self.matrix[a][b].material = 'copperC'
                        self.matrix[a][b].ur = 0.999991
                        self.matrix[a][b].sigma = 5.96*10**6
                    elif b in self.removeLowerCoilIdxs + self.removeUpperCoilIdxs:
                        self.matrix[a][b].material = 'vacuum'
                        self.matrix[a][b].ur = 1.0
                        self.matrix[a][b].sigma = 3*10**(-15)
                    else:
                        self.matrix[a][b].material = ''
                b += 1
            b = 0
            a += 1


class Node(object):

    def __init__(self, iIndex, iXinfo, iZinfo, iDesign):

        self.xIndex = iIndex[0]
        self.yIndex = iIndex[1]

        # Initialize dense meshing near slot and teeth edges in x direction
        self.x = iXinfo[0]
        self.lx = iXinfo[1]

        # Initialize dense meshing near slot and teeth edges in y direction
        self.y = iZinfo[0]
        self.ly = iZinfo[1]

        self.xCenter = self.x + self.lx/2

        self.yCenter = self.y + self.ly/2

        self.ur = 0.0
        self.sigma = 0.0
        self.material = ''
        self.colour = ''

        self.Szy = self.ly*iDesign.D
        self.Sxz = self.lx*iDesign.D

        # Potential
        self.Yk = 0.0

        # Current Density
        self.Jy = 0.0

        # Thrust
        self.Fx, self.Fy, self.F = np.zeros(3, dtype=np.cdouble)

        # Magnetic Field Strength
        self.Hx, self.Hy, self.H = np.zeros(3, dtype=np.cdouble)

        # B field
        self.Bx, self.By, self.B = np.zeros(3, dtype=np.cdouble)

        # Flux
        self.phiXp, self.phiXn, self.phiYp, self.phiYn, self.phiX, self.phiY, self.phi = np.zeros(7, dtype=np.cdouble)
        self.phiError = 0.0

        # Current
        self.Iph = 0.0  # phase

        # MMF
        self.MMF = 0.0  # AmpereTurns

        # Reluctance
        self.Rx, self.Ry, self.R = np.zeros(3, dtype=np.float64)

    '''https://stackoverflow.com/questions/1227121/compare-object-instances-for-equality-by-their-attributes'''
    # This method was created to check if two objects are identical when both use the Node base class.
    # If checking nodeA == nodeB this would return false regardless due to object id but now returns true if equal
    def __eq__(self, otherObject):
        if not isinstance(otherObject, Node):
            # don't attempt to compare against unrelated types
            return NotImplemented
        # If the objects are the same then set the IDs to be equal
        if self.__dict__.items() == otherObject.__dict__.items():
            for attr, val in otherObject.__dict__.items():
                return self.__dict__[attr] == otherObject.__dict__[attr]
        else:
            print('Error - This is saying that the items of one object does not match the other')

    '''https://stackoverflow.com/questions/1227121/compare-object-instances-for-equality-by-their-attributes'''
    # def __hash__(self):
    #     # necessary for instances to behave sanely in dicts and sets.
    #     return hash((self.foo, self.bar))

    def drawNode(self, iNodewidth, gridSpacing, overRideColour, c):

        x_old = self.x*gridSpacing
        x_new = x_old + self.lx*gridSpacing
        y_old = self.y*gridSpacing
        y_new = y_old + self.ly*gridSpacing

        if overRideColour:
            fillColour = overRideColour
        else:
            idx = np.where(matList == self.material)
            if len(idx[0] == 1):
                matIdx = idx[0][0]
                self.colour = matList[matIdx][1]
            else:
                self.colour = 'orange'
            fillColour = self.colour

        c.create_rectangle(x_old, y_old, x_new, y_new, width=iNodewidth, fill=fillColour)

    def methodRebuildA(self, *argv):
        nodeInfo = argv[0]
        attributeNames = list(nodeInfo.keys())
        attributeVals = list(nodeInfo.values())

        for i in range(len(nodeInfo)):
            variableName = attributeNames[i]
            if variableName in self.__dict__:
                # These are all the variables that are complex which JSON can't handle as a dtype
                if type(attributeVals[i]) == list and attributeVals[i][0] == 'plex_Signature':
                    attributeVal = np.cdouble(attributeVals[i][1] + j_plex*attributeVals[i][2])
                else:
                    attributeVal = attributeVals[i]
                self.__dict__[variableName] = attributeVal
        return self


class Region(object):
    def __init__(self, iAn, iBn):
        self.an = iAn
        self.bn = iBn


def reluctance(iNode):
    ResX = iNode.lx/(2*uo*iNode.ur*iNode.Szy)
    ResY = iNode.ly/(2*uo*iNode.ur*iNode.Sxz)
    Res = math.sqrt(ResX ** 2 + ResY ** 2)
    return ResX, ResY, Res


def meshBoundary(iRegionDimensions, iSpacing, iFraction, iNumBoundaries, iMeshDensity):
    meshSize = (iRegionDimensions[0] - iNumBoundaries*iFraction*iSpacing) / (iRegionDimensions[1] - len(iMeshDensity)*iNumBoundaries)
    return meshSize


def checkSpatialMapping(iDesign, iGrid, iPixelDivision, errorList, spatialDomainFlag, iIdxs):

    iIdxLeftAirBuffer, iIdxLeftEndTooth, iIdxSlot, iIdxTooth, iIdxRightEndTooth, iIdxRightAirBuffer = iIdxs

    flag_domainMapping = Error('', False)

    # Check grid spacing to slotpitch
    if round(iDesign.slotpitch - iGrid.Spacing*iPixelDivision, 12) != 0:
        print('flag - iGrid spacing')
        spatialDomainFlag = True
    # Check slotpitch
    if round(iDesign.slotpitch - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxTooth + iGrid.ppTooth].x - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxSlot].x), 12) != 0:
        print('flag - slotpitch')
        spatialDomainFlag = True
    # Check slot width
    if round(iDesign.ws - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxTooth].x - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxSlot].x), 12) != 0:
        print('flag - slots')
        spatialDomainFlag = True
    # Check tooth width
    if round(iDesign.wt - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxTooth + iGrid.ppTooth].x - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxTooth].x), 12) != 0:
        print('flag - teeth')
        spatialDomainFlag = True
    # Check left end tooth
    if round(iDesign.endTeeth - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxSlot].x - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxLeftEndTooth].x), 12) != 0:
        print('flag - left end tooth')
        spatialDomainFlag = True
    # Check right end tooth
    if round(iDesign.endTeeth - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxRightAirBuffer].x - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxRightEndTooth].x), 12) != 0:
        print('flag - right end tooth')
        spatialDomainFlag = True
    # Check left air buffer
    if round(iGrid.AirBuffer - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxLeftEndTooth].x - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxLeftAirBuffer].x), 12) != 0:
        print('flag - left air buffer')
        spatialDomainFlag = True
    # Check right air buffer
    if round(iGrid.AirBuffer - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][-1].x + iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][-1].lx - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][iIdxRightAirBuffer].x), 12) != 0:
        print('flag - right air buffer')
        spatialDomainFlag = True

    # Check Vacuum
    if round(iDesign.vacuumBoundary - (iGrid.matrix[iGrid.ppVacuumLower][0].y - iGrid.matrix[0][0].y), 12) != 0:
        print('flag - vacuum')
        spatialDomainFlag = True
    # Check Yoke
    if round(iDesign.hy - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppYokeheight][0].y - iGrid.matrix[iGrid.ppVacuumLower][0].y), 12) != 0:
        print('flag - yoke')
        spatialDomainFlag = True
    # Check Slot/Tooth Height
    if round(iDesign.hs - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppHeight][0].y - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppHeight-iGrid.ppSlotheight][0].y), 12) != 0:
        print('flag - slot height')
        spatialDomainFlag = True
    # Check Air Gap
    if round(iDesign.g - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppHeight+iGrid.ppAirgap][0].y - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppHeight][0].y), 12) != 0:
        print('flag - air gap')
        spatialDomainFlag = True
    # Check Blade Rotor
    if round(iDesign.dr - (iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppHeight+iGrid.ppAirgap+iGrid.ppBladerotor][0].y - iGrid.matrix[iGrid.ppVacuumLower+iGrid.ppHeight+iGrid.ppAirgap][0].y), 12) != 0:
        print('flag - blade rotor')
        spatialDomainFlag = True

    # This checks that the spatial domain is mapped to the canvas domain
    if spatialDomainFlag:
        flag_domainMapping.description = 'ERROR - The spatial domain does not match with the canvas domain'
        flag_domainMapping.state = True
        errorList.append(flag_domainMapping)

    return None


def GitHub__LIM_Grid(iDesign, iN, iCanvasSpacing, iPixelDivision, errorList, iMeshDensity, iMeshIndexes):

    meshDensity = iMeshDensity
    design = iDesign
    n = iN
    spatialDomainFlag = False
    pixelSpacing = design.slotpitch/iPixelDivision

    grid = Grid(design, pixelSpacing, iCanvasSpacing, meshDensity, n, iMeshIndexes)
    grid.produceGrid(iDesign=design, grid=grid, pixelSpacing=pixelSpacing, iMeshIndexes=iMeshIndexes)

    # Define Indexes
    idxLeftAirBuffer = 0
    idxLeftEndTooth = idxLeftAirBuffer + grid.ppAirBuffer
    idxSlot = idxLeftEndTooth + grid.ppLeftEndTooth
    idxTooth = idxSlot + grid.ppSlot
    idxRightEndTooth = grid.ppL - grid.ppAirBuffer - grid.ppRightEndTooth
    idxRightAirBuffer = idxRightEndTooth + grid.ppRightEndTooth

    # This function is written to catch any errors in the mapping between canvas and space for both x and y coordinates
    checkSpatialMapping(design, grid, iPixelDivision, errorList, spatialDomainFlag, [idxLeftAirBuffer, idxLeftEndTooth, idxSlot, idxTooth, idxRightEndTooth, idxRightAirBuffer])

    windingLayers = 2
    scalingLower, scalingUpper = 0.0, 0.0
    time_plex = cmath.exp(j_plex*2*pi*iDesign.f*iDesign.t)
    i, j = 0, 0
    while i < grid.ppH:
        while j < grid.ppL:
            grid.matrix[i][j].Rx, grid.matrix[i][j].Ry, grid.matrix[i][j].R = reluctance(grid.matrix[i][j])

            if i in grid.lower_slotYIndexes1 and j in grid.coilArray:
                if j in grid.lower_slotsA:
                    angle_plex = cmath.exp(0)
                elif j in grid.lower_slotsB:
                    angle_plex = cmath.exp(-j_plex*pi*2/3)
                elif j in grid.lower_slotsC:
                    angle_plex = cmath.exp(j_plex*pi*2/3)
                else:
                    angle_plex = 0.0

                # Set the scaling factor for MMF in equation 18
                # 2 coils in slot
                if grid.matrix[i][j].material[:-1] == 'copper' and grid.matrix[i-grid.ppSlotheight//2][j].material[:-1] == 'copper':
                    scalingLower = 0.5
                # coil in upper slot only
                elif grid.matrix[i][j].material[:-1] != 'copper' and grid.matrix[i-grid.ppSlotheight//2][j].material[:-1] == 'copper':
                    scalingLower = 0.0
                # coil in lower slot only
                elif grid.matrix[i][j].material[:-1] == 'copper' and grid.matrix[i-grid.ppSlotheight//2][j].material[:-1] != 'copper':
                    scalingLower = 0.5
                # empty slot
                elif grid.matrix[i][j].material == 'vacuum' and grid.matrix[i-grid.ppSlotheight//2][j].material == 'vacuum':
                    scalingLower = 0.0
                else:
                    scalingLower = 0.0

            elif i in grid.upper_slotYIndexes1 and j in grid.coilArray:
                if j in grid.upper_slotsA:
                    angle_plex = cmath.exp(0)
                elif j in grid.upper_slotsB:
                    angle_plex = cmath.exp(-j_plex*pi*2/3)
                elif j in grid.upper_slotsC:
                    angle_plex = cmath.exp(j_plex*pi*2/3)
                else:
                    angle_plex = 0.0

                # Set the scaling factor for MMF in equation 18
                # 2 coils in slot
                if grid.matrix[i][j].material[:-1] == 'copper' and grid.matrix[i+grid.ppSlotheight//2][j].material[:-1] == 'copper':
                    scalingUpper = 1.5
                # coil in lower slot only
                elif grid.matrix[i][j].material[:-1] != 'copper' and grid.matrix[i+grid.ppSlotheight//2][j].material[:-1] == 'copper':
                    scalingUpper = 1.0
                # coil in upper slot only
                elif grid.matrix[i][j].material[:-1] == 'copper' and grid.matrix[i+grid.ppSlotheight//2][j].material[:-1] != 'copper':
                    scalingUpper = 0.5
                # empty slot
                elif grid.matrix[i][j].material == 'vacuum' and grid.matrix[i+grid.ppSlotheight//2][j].material == 'vacuum':
                    scalingUpper = 0.0
                else:
                    scalingUpper = 0.0
            else:
                angle_plex = 0.0
                scalingLower = 0.0
                scalingUpper = 0.0

            # turnAreaRatio is the fractional area of MEC nodes per coil
            turnAreaRatio = grid.matrix[i][j].lx * grid.matrix[i][j].ly / (windingLayers * iDesign.hs * iDesign.ws)
            grid.matrix[i][j].Iph = (iDesign.Iin*math.sqrt(2))*angle_plex*time_plex

            # Lower slots only
            if i in grid.lower_slotYIndexes1 and j in grid.coilArray:
                if j in grid.inLower_slotsA or j in grid.inLower_slotsB or j in grid.inLower_slotsC:
                    inOutCoeffMMF = -1
                elif j in grid.outLower_slotsA or j in grid.outLower_slotsB or j in grid.outLower_slotsC:
                    inOutCoeffMMF = 1
                else:
                    if j not in grid.removeLowerCoilIdxs:
                        print('Shouldnt happen')
                    inOutCoeffMMF = 0

                grid.matrix[i][j].MMF = inOutCoeffMMF * scalingLower * iDesign.N * grid.matrix[i][j].Iph / (2 * turnAreaRatio)

            # Upper slots only
            elif i in grid.upper_slotYIndexes1 and j in grid.coilArray:
                if j in grid.inUpper_slotsA or j in grid.inUpper_slotsB or j in grid.inUpper_slotsC:
                    inOutCoeffMMF = -1
                elif j in grid.outUpper_slotsA or j in grid.outUpper_slotsB or j in grid.outUpper_slotsC:
                    inOutCoeffMMF = 1
                else:
                    if j not in grid.removeUpperCoilIdxs:
                        print('Shouldnt happen')
                    inOutCoeffMMF = 0

                grid.matrix[i][j].MMF = inOutCoeffMMF * scalingUpper * iDesign.N * grid.matrix[i][j].Iph / (2 * turnAreaRatio)

            # Stator teeth
            else:
                grid.matrix[i][j].MMF = 0.0

            j += 1
        j = 0
        i += 1

    return grid, errorList
