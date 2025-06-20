# add testing subclasses here

from path_utils import path_to_fates_root, prepend_to_python_path
import os
_FUNCTIONAL_TESTING_DIR = os.path.join(path_to_fates_root(), "testing")
prepend_to_python_path(_FUNCTIONAL_TESTING_DIR)

from functional_class import FunctionalTest
from functional_testing.allometry.allometry_test import AllometryTest
from functional_testing.math_utils.math_utils_test import QuadraticTest
from functional_testing.fire.fuel.fuel_test import FuelTest 
from functional_testing.fire.ros.ros_test import ROSTest  
from functional_testing.patch.patch_test import PatchTest
