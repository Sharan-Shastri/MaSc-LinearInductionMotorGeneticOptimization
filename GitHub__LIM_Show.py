from GitHub__LIM_Grid import *
import tkinter as tk
import matplotlib as mpl


# noinspection PyUnresolvedReferences
def colorFader(c1, c2, mix=0):  # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)

    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)


def myColourNumber(iFieldScale, iVal):
    if iVal in [np.inf, -np.inf]:
        result = iVal
    else:
        result = min(iFieldScale, key=lambda x: abs(x - iVal))

    return result


def determineColour(iGrid, iGridInfo, iI, iJ, iField):

    iFieldType, iFieldsScale, iStoColours, iPosInf, iNegInf = iField

    # I just added this, Might need to add it elsewhere. Make sure its correct
    if type(iGrid[iI, iJ].__dict__[iFieldType]) in iGridInfo['complexTypeList']:
        myNumber = myColourNumber(iFieldsScale, iGrid[iI, iJ].__dict__[iFieldType].real)
    else:
        myNumber = myColourNumber(iFieldsScale, iGrid[iI, iJ].__dict__[iFieldType])

    if myNumber == np.inf:
        oOverRideColour = iPosInf
    elif myNumber == -np.inf:
        oOverRideColour = iNegInf
    else:
        # noinspection PyUnboundLocalVariable
        colorScaleIndex = np.where(iFieldsScale == myNumber)
        if colorScaleIndex[0][0] == len(iStoColours):
            print('OH MY LAWD ITS A FIRE')
            colorScaleIndex[0][0] = len(iStoColours) - 1
        oOverRideColour = iStoColours[colorScaleIndex[0][0]]

    return oOverRideColour


def minMaxField(iInfo, iGrid, iAttName, iFiltered, iShowFilter):

    iFilteredRows, iFilteredRowCols = iFiltered

    if iShowFilter:
        tFiltered = [[iGrid[x.yIndex, x.xIndex] for x in y if x.yIndex in iFilteredRows or (x.yIndex, x.xIndex) in iFilteredRowCols] for y in iGrid]

    else:
        tFiltered = iGrid

    tFilteredNoEmpties = list(filter(lambda x: True if x != [] else False, tFiltered))

    maxAttr = max(map(max, [[x.__dict__[iAttName].real if str(type(x.__dict__[iAttName])).split("'")[1] in iInfo['complexTypeList'] else x.__dict__[iAttName] for x in y] for y in tFilteredNoEmpties]))
    minAttr = min(map(min, [[x.__dict__[iAttName].real if str(type(x.__dict__[iAttName])).split("'")[1] in iInfo['complexTypeList'] else x.__dict__[iAttName] for x in y] for y in tFilteredNoEmpties]))
    maxScale = maxAttr
    minScale = minAttr

    return minScale, maxScale


def combineFilterList(i_pp, iUnfilteredRows, iUnfilteredRowCols):

    i_ppH, i_ppL = i_pp
    oFilteredRows = []
    oFilteredRowCols = []

    # List of row indexes for keeping rows from 0 to ppH
    for a in np.arange(len(iUnfilteredRows)):
        oFilteredRows += iUnfilteredRows[a]

    # List of (row, col) indexes for keeping nodes in mesh
    for a in np.arange(len(iUnfilteredRowCols)):
        for y in np.arange(i_ppH):
            for x in np.arange(i_ppL):
                if y in iUnfilteredRowCols[a][0] and x in iUnfilteredRowCols[a][1]:
                    oFilteredRowCols += [(y, x)]

    return oFilteredRows, oFilteredRowCols


def GitHub__LIM_Show(iGridInfo, iGridMatrix, iFieldType, iShowGrid, iShowFields, iShowFilter, iNumColours, errorList, iDims, iFieldsList):

    grid = iGridMatrix
    gridInfo = iGridInfo
    fieldType = iFieldType
    showGrid = iShowGrid
    showFields = iShowFields
    showFilter = iShowFilter
    dims = iDims
    fieldsList = iFieldsList

    flag_wrongField, flag_emptyField = Error('', False), Error('', False)

    # Create the grid canvas to display the grid mesh
    if showGrid:
        rootGrid = tk.Tk()
        cGrid: tk.Canvas = tk.Canvas(rootGrid, height=dims[0], width=dims[1], bg='gray30')
        cGrid.pack()
        i, j = 0, 0
        while i < grid.shape[0]:
            while j < grid.shape[1]:
                grid[i, j].drawNode(iNodewidth=1, gridSpacing=gridInfo['Cspacing'], overRideColour=False, c=cGrid)
                j += 1
            j = 0
            i += 1
        rootGrid.mainloop()

    if fieldType not in fieldsList:
        flag_wrongField.description = 'WARNING - Field Analysis Error. Selected Field Type Not In Fields List %s' % fieldsList
        flag_wrongField.state = True
        errorList.append(flag_wrongField)

    # Color nodes based on node values
    elif showFields or showFilter:

        rootFields = tk.Tk()
        cFields: tk.Canvas = tk.Canvas(rootFields, height=dims[0], width=dims[1], bg='gray30')
        cFields.pack()

        c1 = '#FFF888'  # Yellow Positive Limit
        c2 = '#700000'  # Dark Red Negative Limit
        cPosInf = '#00FFFF'  # BabyBlue Positive Infinity
        cNegInf = '#9EFE4C'  # Green Negative Infinity
        stoColours = ["empty"]*(iNumColours + 1)
        for x in range(iNumColours + 1):
            stoColours[x] = colorFader(c1, c2, x / iNumColours)

        # Max and Min values for normalizing the color scales for field analysis
        keepRows = [[]]*1
        keepRowColsUnfiltered = [[]]*2

        # [row]
        # Rule 1
        keepRows[0] = gridInfo['airgapYIndexes']

        # [row, col] - Make sure to put the rules in order of ascending rows or the list wont be sorted (shouldnt matter)
        # Rule 1
        keepRowColsUnfiltered[0] = [gridInfo['yokeYIndexes'], gridInfo['toothArray'] + gridInfo['coilArray']]
        # Rule 2
        keepRowColsUnfiltered[1] = [gridInfo['lower_slotYIndexes1'] + gridInfo['upper_slotYIndexes1'], gridInfo['toothArray']]

        filteredRows, filteredRowCols = combineFilterList([gridInfo['ppH'], gridInfo['ppL']], keepRows, keepRowColsUnfiltered)

        minScale, maxScale = minMaxField(gridInfo, grid, fieldType, [filteredRows, filteredRowCols], iShowFilter=showFilter)
        normScale = (maxScale - minScale) / (iNumColours - 1)

        if [minScale, maxScale, normScale] == [0, 0, 0]:
            flag_emptyField.description = 'WARNING - Field Analysis Error. All values are zero! Type: ' + fieldType
            flag_emptyField.state = True
            errorList.append(flag_emptyField)

        # Create fields canvas to display the selected field result on the mesh
        else:
            if minScale < 0:
                fieldsScale = np.arange(minScale, maxScale - normScale, normScale)
            else:
                fieldsScale = np.arange(minScale, maxScale + normScale, normScale)
            colorScaleIndex = np.where(fieldsScale == fieldsScale[0])

            if minScale < 0:
                cFields.create_text(400, 1000, font="Purisa", text="Debug (Max, Min): (%s, %s) colour: (%s, %s) Type: %s" % (maxScale, minScale, stoColours[-1], stoColours[colorScaleIndex[0][0]], fieldType))
                print("Debug (Max, Min): (%s, %s) colour: (%s, %s) Type: %s" % (maxScale, minScale, stoColours[-1], stoColours[colorScaleIndex[0][0]], fieldType))
            else:
                cFields.create_text(400, 1000, font="Purisa", text="Debug (Max, Min): (%s, %s) colour: (%s, %s) Type: %s" % (maxScale, minScale, stoColours[-1], stoColours[colorScaleIndex[0][0]], fieldType))
                print("Debug (Max, Min): (%s, %s) colour: (%s, %s) Type: %s" % (maxScale, minScale, stoColours[-1], stoColours[colorScaleIndex[0][0]], fieldType))

            # Assigns a colour to a node based on its relative position in the range of values and the range of available colours
            i, j, k = 0, 0, 0
            while i < grid.shape[0]:
                while j < grid.shape[1]:
                    if showFilter:
                        if i in filteredRows:
                            overRideColour = determineColour(grid, gridInfo, i, j, [fieldType, fieldsScale, stoColours, cPosInf, cNegInf])

                        elif (i, j) in filteredRowCols:
                            overRideColour = determineColour(grid, gridInfo, i, j, [fieldType, fieldsScale, stoColours, cPosInf, cNegInf])

                        else:
                            overRideColour = '#000000'
                    else:
                        overRideColour = determineColour(grid, gridInfo, i, j, [fieldType, fieldsScale, stoColours, cPosInf, cNegInf])

                    grid[i, j].drawNode(iNodewidth=1, gridSpacing=gridInfo['Cspacing'], overRideColour=overRideColour, c=cFields)
                    j += 1
                    k += 1
                j = 0
                i += 1

        rootFields.mainloop()

    return fieldsList, errorList
