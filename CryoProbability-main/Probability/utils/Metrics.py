"""
Metrics used when calculating the coverage of the iceball around the tumor
    - Percentage Coverage
    - Probability
    - Dice Coefficient Score
    - Target Overlap
    - Positive Predictive Value
    - Time
"""
class quad:

    """def createPerturbedCenters(self):
        num_nodes = 5
        perturbed_x, weights = self.gaussian_quadrature_nodes(num_nodes, center[0], SDinplane)
        perturbed_y, weights = self.gaussian_quadrature_nodes(num_nodes, center[1], SDinplane)
        perturbed_z, weights = self.gaussian_quadrature_nodes(num_nodes, center[2], SDdepth)"""

    def runQuadrature(self, center, inputVolume, SDinplane, SDdepth, nOfPoints):
        num_nodes = 5
        coordinates = []
        coord = numpy.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
        for i in range(nOfPoints):
            perturbed_x, weights = self.gaussian_quadrature_nodes(num_nodes, center[i][0], SDinplane)
            perturbed_y, weights = self.gaussian_quadrature_nodes(num_nodes, center[i][1], SDinplane)
            perturbed_z, weights = self.gaussian_quadrature_nodes(num_nodes, center[i][2], SDdepth)
            coordinate = self.create_coordinates_new(perturbed_x, perturbed_y, perturbed_z)
            for j in range(num_nodes):
                coord[i][0 + j * 3] = coordinate[j][0]
                coord[i][1 + j * 3] = coordinate[j][1]
                coord[i][2 + j * 3] = coordinate[j][2]
        percentage, prob, dice_coefficient_score, target_overlap, ppv_score = self.perturbed_tumor_percentage(
            inputVolume, coord, [ICESEED_R, ICESEED_A, ICESEED_S], weights, nOfPoints
        )
        return percentage, prob, dice_coefficient_score, target_overlap, ppv_score

    def gaussian_quadrature_nodes_weights(self, n):
        # Generate nodes and weights for Legendre-Gauss quadrature
        nodes, weights = numpy.polynomial.legendre.leggauss(n)

        # nodes, weights = numpy.polynomial.hermite.hermgauss(n)
        # Map nodes to the interval [-1, 1]
        nodes = nodes * 2 - 1
        return nodes, weights

    def create_coordinates_new(self, x_list, y_list, z_list):
        # Combine all elements from the lists without repetition
        coordinates = [(x_list[0], y_list[0], z_list[0]), (x_list[1], y_list[1], z_list[1]),
                       (x_list[2], y_list[2], z_list[2]), (x_list[3], y_list[3], z_list[3]),
                       (x_list[4], y_list[4], z_list[4])]
        return coordinates

    def create_coordinates(self, x_list, y_list, z_list):
        # Combine all elements from the lists without repetition
        coordinates = list(product(x_list, y_list, z_list))
        return coordinates

    def perturbed_tumor_percentage(self, nrrd, coordinates, radii, weights, nOfProbes):
        tumor_percentages = []
        prob = []
        dice_coefficient_scores = []
        target_overlaps = []
        ppv_scores = []

        for j in range(5):  # number of nodes - change
            center = [[coordinates[0][0 + j * 3], coordinates[0][1 + j * 3], coordinates[0][2 + j * 3]],
                      [coordinates[1][0 + j * 3], coordinates[1][1 + j * 3], coordinates[1][2 + j * 3]],
                      [coordinates[2][0 + j * 3], coordinates[2][1 + j * 3], coordinates[2][2 + j * 3]]]

            percentage = self.iceball_coverage_1(nrrd, center, nOfProbes)
            dice_coefficient_score = self.dice_coefficient(nrrd, center, nOfProbes)
            target_overlap = self.target_overlap(nrrd, center, nOfProbes)
            ppv_score = self.positive_predictive_value(nrrd, center, nOfProbes)

            tumor_percentages.append(percentage * weights[j])
            dice_coefficient_scores.append(dice_coefficient_score * weights[j])
            target_overlaps.append(target_overlap * weights[j])
            ppv_scores.append(ppv_score * weights[j])

            if percentage >= 0.99:
                prob.append(weights[j])
            else:
                prob.append(0.0)

        return(
            numpy.sum(tumor_percentages) / numpy.sum(weights),
            numpy.sum(prob) / numpy.sum(weights),
            numpy.sum(dice_coefficient_scores) / numpy.sum(weights),
            numpy.sum(target_overlaps) / numpy.sum(weights),
            numpy.sum(ppv_scores) / numpy.sum(weights)
        )

    def gaussian_quadrature_nodes(self, num_nodes, center, sd):
        nodes, weights = self.gaussian_quadrature_nodes_weights(num_nodes)
        # Perturb each node value by the specified standard deviation
        perturbed_nodes = center + nodes * sd
        # Return the perturbed node values
        return perturbed_nodes, weights


    def getIceballImageData(self, center, Spacing, size_image):
        sphere1 = vtk.vtkImageEllipsoidSource()
        sphere1.SetOutputScalarTypeToShort()
        sphere1.SetCenter([center[0], center[1], center[2] - 3])
        sphere1.SetRadius(ICESEED_R / Spacing[0], ICESEED_A / Spacing[1], ICESEED_S / Spacing[2])
        sphere1.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
        sphere1.Update()
        return sphere1.GetOutput()


    def andLogic(self, image, sp):
        logic = vtk.vtkImageLogic()
        logic.SetInput1Data(image)
        logic.SetInput2Data(sp)
        logic.SetOperationToAnd()
        logic.SetOutputTrueValue(1)
        logic.Update()
        return logic.GetOutput()


    def orLogic(self, image, sp):
        logic = vtk.vtkImageLogic()
        logic.SetInput1Data(image)
        logic.SetInput2Data(sp)
        logic.SetOperationToOr()
        logic.SetOutputTrueValue(1)
        logic.Update()
        return logic.GetOutput()


    def addIceballToScene(self, image, IjkToRasMatrix, data):
        try:
            outputLabelmapVolumeNode = slicer.util.getNode('iceball')
        except:
            outputLabelmapVolumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode")
            outputLabelmapVolumeNode.SetName("iceball")

        outputLabelmapVolumeNode.SetOrigin(image.GetOrigin())
        outputLabelmapVolumeNode.SetSpacing(image.GetSpacing())
        outputLabelmapVolumeNode.SetIJKToRASMatrix(IjkToRasMatrix)
        outputLabelmapVolumeNode.SetAndObserveImageData(data)
        outputLabelmapVolumeNode.CreateDefaultDisplayNodes()
        slicer.util.setSliceViewerLayers(outputLabelmapVolumeNode, fit=True)


    def iceball_coverage_1(self, reader, center, nOfProbes):
        RASToIJKMatrix = vtk.vtkMatrix4x4()
        reader.GetRASToIJKMatrix(RASToIJKMatrix)
        IjkToRasMatrix = vtk.vtkMatrix4x4()
        reader.GetIJKToRASMatrix(IjkToRasMatrix)

        image = reader.GetImageData()
        image.AllocateScalars(4, 1)
        size_image = image.GetDimensions()
        Spacing = reader.GetSpacing()

        position = [center[0][0], center[0][1], center[0][2], 1]
        center_probe = RASToIJKMatrix.MultiplyPoint(position)
        sp = self.getIceballImageData(center_probe, Spacing, size_image)

        for i in range(nOfProbes):
            position = [center[i][0], center[i][1], center[i][2], 1]
            center_probe = RASToIJKMatrix.MultiplyPoint(position)
            sp2 = self.getIceballImageData(center_probe, Spacing, size_image)
            sp = self.orLogic(sp, sp2)

        self.addIceballToScene(reader, IjkToRasMatrix, sp)

        fusedImage = self.andLogic(image, sp)

        # Extract scalar data and calculate iceball coverage
        image_data = vtk_to_numpy(image.GetPointData().GetScalars())
        logic_data = vtk_to_numpy(fusedImage.GetPointData().GetScalars())

        image_coverage = numpy.sum(image_data)  # Total coverage of the original image
        iceball_coverage = numpy.sum(logic_data)  # Coverage of the iceball within the image

        return iceball_coverage / image_coverage

    # 2(S∩Σ) /(S+Σ)
    def dice_coefficient(self, reader, center, nOfProbes):
        RASToIJKMatrix = vtk.vtkMatrix4x4()
        reader.GetRASToIJKMatrix(RASToIJKMatrix)
        IjkToRasMatrix = vtk.vtkMatrix4x4()
        reader.GetIJKToRASMatrix(IjkToRasMatrix)

        image = reader.GetImageData()
        image.AllocateScalars(4, 1)
        size_image = image.GetDimensions()
        Spacing = reader.GetSpacing()

        sp = self.getIceballImageData(center[0], Spacing, size_image)
        for i in range(nOfProbes):
            position = [center[i][0], center[i][1], center[i][2], 1]
            center_probe = RASToIJKMatrix.MultiplyPoint(position)
            sp2 = self.getIceballImageData(center_probe, Spacing, size_image)
            sp = self.orLogic(sp, sp2)

        self.addIceballToScene(reader, IjkToRasMatrix, sp)

        fusedImage = self.andLogic(image, sp)

        # Calculate Dice Coefficient
        intersection = vtk_to_numpy(fusedImage.GetPointData().GetScalars())
        dice = 2 * numpy.sum(intersection) / (numpy.sum(vtk_to_numpy(sp.GetPointData().GetScalars())) + numpy.sum(
            vtk_to_numpy(image.GetPointData().GetScalars())))
        return dice

        # (|S∩Σ|)/(|S|)
    def target_overlap(self, reader, center, nOfProbes):
        RASToIJKMatrix = vtk.vtkMatrix4x4()
        reader.GetRASToIJKMatrix(RASToIJKMatrix)
        IjkToRasMatrix = vtk.vtkMatrix4x4()
        reader.GetIJKToRASMatrix(IjkToRasMatrix)

        image = reader.GetImageData()
        image.AllocateScalars(4, 1)
        size_image = image.GetDimensions()
        Spacing = reader.GetSpacing()

        sp = self.getIceballImageData(center[0], Spacing, size_image)
        for i in range(nOfProbes):
            position = [center[i][0], center[i][1], center[i][2], 1]
            center_probe = RASToIJKMatrix.MultiplyPoint(position)
            sp2 = self.getIceballImageData(center_probe, Spacing, size_image)
            sp = self.orLogic(sp, sp2)

        self.addIceballToScene(reader, IjkToRasMatrix, sp)

        fusedImage = self.andLogic(image, sp)

        # Calculate Target Overlap
        intersection = vtk_to_numpy(fusedImage.GetPointData().GetScalars())
        overlap = numpy.sum(intersection) / numpy.sum(vtk_to_numpy(image.GetPointData().GetScalars()))
        return overlap

        # (|S∩Σ|)/(|Σ|
    def positive_predictive_value(self, reader, center, nOfProbes):
        RASToIJKMatrix = vtk.vtkMatrix4x4()
        reader.GetRASToIJKMatrix(RASToIJKMatrix)
        IjkToRasMatrix = vtk.vtkMatrix4x4()
        reader.GetIJKToRASMatrix(IjkToRasMatrix)

        image = reader.GetImageData()
        image.AllocateScalars(4, 1)
        size_image = image.GetDimensions()
        Spacing = reader.GetSpacing()

        sp = self.getIceballImageData(center[0], Spacing, size_image)
        for i in range(nOfProbes):
            position = [center[i][0], center[i][1], center[i][2], 1]
            center_probe = RASToIJKMatrix.MultiplyPoint(position)
            sp2 = self.getIceballImageData(center_probe, Spacing, size_image)
            sp = self.orLogic(sp, sp2)

        self.addIceballToScene(reader, IjkToRasMatrix, sp)

        fusedImage = self.andLogic(image, sp)

        # Calculate Positive Predictive Value
        intersection = vtk_to_numpy(fusedImage.GetPointData().GetScalars())
        ppv = numpy.sum(intersection) / numpy.sum(vtk_to_numpy(sp.GetPointData().GetScalars()))
        return ppv