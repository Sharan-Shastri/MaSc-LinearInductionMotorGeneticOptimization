from GitHub__LIM_Show import *
import json
import inquirer


def userInput(iFieldsList):

    # Resolution
    qReso = [
        inquirer.List('resolution',
                      message='Enter Monitor Resolution',
                      choices=[[1280, 1024], [1366, 768], [1600, 900], [1920, 1080], [1920, 1200], [2560, 1440], [3440, 1440], [3840, 2160]],
                      ),
    ]
    answers = inquirer.prompt(qReso)
    resolution = answers['resolution']

    # Show Grid
    confirmGrid = {
        inquirer.Confirm('grid',
                         message='Do You Want To Show The Grid?',
                         default=True),
    }
    confirmation = inquirer.prompt(confirmGrid)
    bGrid = confirmation['grid']

    # Show Fields
    confirmFields = {
        inquirer.Confirm('fields',
                         message='Do You Want To Show The Fields Plot?',
                         default=True),
    }
    confirmation = inquirer.prompt(confirmFields)
    bFields = confirmation['fields']

    # Show Filter
    confirmFilter = {
        inquirer.Confirm('filter',
                         message='Do You Want To Show The Fields Filer?',
                         default=True),
    }
    confirmation = inquirer.prompt(confirmFilter)
    bFilter = confirmation['filter']

    # Fields Attribute
    qAttr = [
        inquirer.List('attribute',
                      message='Which Attribute Do You Want To Show In Fields %s?' % iFieldsList,
                      choices=iFieldsList,
                      ),
    ]
    answers = inquirer.prompt(qAttr)
    attribute = answers['attribute']

    # Show Extended Options
    confirmOptions = {
        inquirer.Confirm('options',
                         message='Do You Want To Extended Options?',
                         default=True),
    }

    confirmation = inquirer.prompt(confirmOptions)
    bExtendedOptions = confirmation['options']

    if bExtendedOptions:
        # Mesh Density
        qDensity = [
            inquirer.List('density',
                          message='Enter Mesh Density',
                          choices=[2, 5, 10],
                          ),
        ]
        answers = inquirer.prompt(qDensity)
        density = answers['density']

        # Slots
        qSlots = [
            inquirer.List('slots',
                          message='Enter Number Of Slots',
                          choices=[16, 18, 25, 30],
                          ),
        ]
        answers = inquirer.prompt(qSlots)
        slots = answers['slots']

        # Poles
        qPoles = [
            inquirer.List('poles',
                          message='Enter Number Of Poles',
                          choices=[2, 4, 6, 8],
                          ),
        ]
        answers = inquirer.prompt(qPoles)
        poles = answers['poles']

        dictOptions = {'density': density, 'slots': slots, 'poles': poles}

    else:
        dictOptions = {'density': False, 'slots': False, 'poles': False}

    # Execute Program
    confirmExecute = {
        inquirer.Confirm('execute',
                         message='Are You Ready To Run The Program?',
                         default=True),
    }

    confirmation = inquirer.prompt(confirmExecute)
    bExecute = confirmation['execute']

    returnInputs = {'resolution': resolution, 'grid': bGrid, 'fields': bFields, 'filter': bFilter,
                    'attribute': attribute, 'options': bExtendedOptions, 'optionsDict': dictOptions,
                    'execute': bExecute}

    return returnInputs


# Break apart the grid.matrix array
def destruct(iGrid):
    matrix = iGrid.matrix
    dictList = []
    tTypeList = []

    for row in matrix:
        for col in row:
            nodeList = []
            for attr, val in col.__dict__.items():
                # Generate a list of data types before destructing the matrix
                if str(type(val)).split("'")[1] not in tTypeList:
                    tStrType = str(type(val)).split("'")
                    tTypeList += [tStrType[1]]
                # JSON dump cannot handle numpy arrays
                if type(val) == np.ndarray:
                    val = list(val)
                # Json dump and load cannot handle complex objects and must be an accepted dtype such as list
                elif type(val) in [complex, np.complex128]:
                    # 'plex_Signature' is used to identify if a list is a destructed complex number or not
                    val = ['plex_Signature', val.real, val.imag]
                else:
                    pass

                nodeList.append((attr, val))
            x = dict(nodeList)
            dictList.append(x)

    nonComplexList = ['int', 'np.float64', 'float', 'str']
    tComplexList = [i for i in tTypeList if i not in nonComplexList]
    iGrid.typeList = tTypeList
    iGrid.complexTypeList = tComplexList

    # InfoRes includes all the grid info (excluding the matrix) that is serializable for a JSON object
    infoRes = {attr: val for attr, val in iGrid.__dict__.items() if attr != 'matrix' and type(val) != np.ndarray and str(type(val)).split("'")[1] not in iGrid.complexTypeList}
    # MatrixRes is the destructed matrix array from the grid
    matrixRes = {i: dictList[i] for i in range(len(dictList))}
    dictRes = {'info': infoRes, 'matrix': matrixRes}

    return dictRes


# Rebuild the grid.matrix array
def construct(iDict, iArrayShape, iDesign):
    info = iDict['info']
    matrix = iDict['matrix']
    lenKeys = len(dict.keys(matrix))
    B = np.array([type('', (object,), {}) for x in np.arange(lenKeys)])

    for key in dict.keys(matrix):
        nodeInfo = matrix[key]
        emptyNode = Node([0, 0], [0, 0], [0, 0], iDesign)
        rebuiltNode = emptyNode.methodRebuildA(nodeInfo)
        B[int(key)] = rebuiltNode
    B = B.reshape(iArrayShape[0], iArrayShape[1])
    return info, B


def jsonStoreSolution(iGrid):
    destructedMatA = destruct(iGrid)
    # Data written to file
    with open("StoredSolutionData.json", 'w') as StoredSolutionData:
        json.dump(destructedMatA, StoredSolutionData)


def jsonRestoreSolution(iGrid, iMotor):
    # Data read from file
    with open("StoredSolutionData.json") as StoredSolutionData:
        dictionary = json.load(StoredSolutionData)

    gridInfo, rebuilt__Grid_matrix = construct(iDict=dictionary, iArrayShape=iGrid.matrix.shape, iDesign=iMotor)
    checkIdenticalLists = np.array([[iGrid.matrix[y, x] == rebuilt__Grid_matrix[y, x] for x in np.arange(iGrid.ppL)] for y in np.arange(iGrid.ppH)])
    return gridInfo, rebuilt__Grid_matrix, checkIdenticalLists


def main(iFieldsList, iInputs):
    with timing():

        flag_gridMatrixJSON, flag_MeshDensityDiscrepancy = Error('', False), Error('', False)

        dims = iInputs['resolution']
        dims.reverse()
        lowDiscrete = 50
        n = np.arange(-lowDiscrete, lowDiscrete + 1, dtype=np.int32)

        if iInputs['options']:
            pixelDivisions = iInputs['optionsDict']['density']
        else:
            pixelDivisions = 5

        if iInputs['options']:
            slots = iInputs['optionsDict']['slots']
        else:
            slots = 16
        if iInputs['options']:
            poles = iInputs['optionsDict']['poles']
        else:
            poles = 6

        wt, ws = 6 / 1000, 10 / 1000
        slotpitch = wt + ws
        endTeeth = 2 * (1.5 * wt)
        length = ((slots - 1) * slotpitch + ws) + endTeeth

        errorList = []
        tempMotor = LimMotor(iSlots=slots, iPoles=poles, iL=length)

        # This value defines how small the mesh is at [border, border+1]. Ex) [4, 2] means that the mesh at the border is 1/4 the mesh far away from the border
        meshDensity = np.array([4, 2])
        # Change the mesh density at boundaries. A 1 indicates denser mesh at a size of len(meshDensity)
        # [LeftAirBuffer], [LeftEndTooth], [Slots], [FullTeeth], [LastSlot], [RightEndTooth], [RightAirBuffer]
        xMeshIndexes = [[0, 0]] + [[0, 0]] + [[0, 0], [0, 0]] * (tempMotor.slots - 1) + [[0, 0]] + [[0, 0]] + [[0, 0]]
        # [LowerVac], [Yoke], [LowerSlots], [UpperSlots], [Airgap], [BladeRotor], [BackIron], [UpperVac]
        yMeshIndexes = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]

        if xMeshIndexes[2] != xMeshIndexes[-3]:
            flag_MeshDensityDiscrepancy.description = 'ERROR - The last slot has a different mesh density than all other slots'
            flag_MeshDensityDiscrepancy.state = True
            errorList.append(flag_MeshDensityDiscrepancy)

        tempGrid, errorList = GitHub__LIM_Grid(iDesign=tempMotor, iN=n, iCanvasSpacing=100 / tempMotor.H,
                                               iPixelDivision=pixelDivisions, errorList=errorList,
                                               iMeshDensity=meshDensity,
                                               iMeshIndexes=[xMeshIndexes, yMeshIndexes])

        # Pass contents to json for file storage/manipulation
        jsonStoreSolution(tempGrid)

        # Return from json to display results
        grid_info, grid_matrix, boolIdenticalLists = jsonRestoreSolution(tempGrid, tempMotor)

        if np.all(np.all(boolIdenticalLists, axis=1)):
            # iDims (height x width): BenQ = 1440 x 2560, ViewSonic = 1080 x 1920
            iFieldsList, errorList = GitHub__LIM_Show(grid_info, grid_matrix, iFieldType=iInputs['attribute'], iShowGrid=iInputs['grid'], iShowFields=iInputs['fields'],
                                         iShowFilter=iInputs['filter'], iNumColours=350, errorList=errorList, iDims=dims, iFieldsList=iFieldsList)

        else:
            flag_gridMatrixJSON.description = 'ERROR - The JSON object matrix does not match the original matrix'
            flag_gridMatrixJSON.state = True
            errorList.append(flag_gridMatrixJSON)

        if errorList:
            print()
            print('Below are a list of warnings and errors:')
            for error in errorList:
                print(error.description)


if __name__ == '__main__':
    '''User Input'''

    fieldsList = ['yCenter', 'xCenter', 'Rx', 'Ry', 'R', 'MMF', 'Yk', 'Iph', 'Szy', 'Sxz', 'Bx', 'By', 'B', 'phiXp', 'phiXn', 'phiYp', 'phiYn', 'phiX', 'phiY', 'phi', 'phiError', 'Fx', 'Fy', 'F']

    inputs = userInput(fieldsList)

    if inputs['execute']:
        main(fieldsList, inputs)
    else:
        print('Please Rethink Your Choices :(')

