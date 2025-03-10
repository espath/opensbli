from sympy import Float

class SimulationDataType(object):
    @staticmethod
    def set_datatype(types):
        if types == Float:
            raise ValueError("Incorrect Simulation DataType. Please use FloatC, not Float.")
        else:
            SimulationDataType.dtype = types

    @staticmethod
    def dtype():
        return SimulationDataType.dtype

    @staticmethod
    def opsc():
        return SimulationDataType.dtype.opsc()


class DataType(object):
    pass


class Double(DataType):

    @staticmethod
    def opsc():
        return "double"


class FloatC(DataType):
    @staticmethod
    def opsc():
        return "float"

class Half(DataType):
    @staticmethod
    def opsc():
        return "half"


class UserDefined(DataType):
    """ User defined datatype this is either float or double depending on input"""

    def __init__(self):
        return


class Int(DataType):
    @staticmethod
    def opsc():
        return "int"
