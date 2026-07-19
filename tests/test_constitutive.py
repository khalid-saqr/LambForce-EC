from lambforce_ec.constitutive import sls_complex_modulus
from lambforce_ec.models import SLSMaterial


def test_sls_is_passive_and_has_correct_limits():
    material = SLSMaterial(1000.0, 2.0, 0.1)
    assert sls_complex_modulus(material, 0.0) == 1000.0
    value = sls_complex_modulus(material, 10.0)
    assert value.real > 1000.0
    assert value.imag > 0.0
    high = sls_complex_modulus(material, 1e12)
    assert abs(high.real / 2000.0 - 1) < 1e-8
