import numpy

class PlotData2D(object):
    def __init__(self, _xAxis, _yAxis, _data,
                 _xAxisName, _yAxisName, _dataName) :
        self.xAxis = _xAxis
        self.yAxis = _yAxis
        self.data = _data
        self.__xAxisName = _xAxisName
        self.__yAxisName = _yAxisName
        self.__dataName = _dataName

    def xAxisName(self) :
        return self.__xAxisName

    def yAxisName(self) :
        return self.__yAxisName

    def dataName(self) :
        return self.__dataName


class PlotData1D(object):
    def __init__(self, _xAxis, _data,
                 _xAxisName, _yAxisName, _dataName) :
        self.xAxis = _xAxis
        self.data = _data
        self.__xAxisName = _xAxisName
        self.__yAxisName = _yAxisName
        self.__dataName = _dataName

    def xAxisName(self) :
        return self.__xAxisName

    def yAxisName(self) :
        return self.__yAxisName

    def dataName(self) :
        return self.__dataName
        
