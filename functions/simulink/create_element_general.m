function element = create_element_general(matrix, name)
    element = Simulink.BusElement;
    element.Name = name;
    element.Dimensions = size(matrix);
    element.DimensionsMode = 'Fixed';
    element.DataType = 'double';
    element.SampleTime = -1;
    element.Complexity = 'real';
end
