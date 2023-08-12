import logging
import kpp_to_micm

logger = logging.getLogger(__name__)

def test_parse_kpp_arrhenius():
    """
    examples
    ARR_ab(1.0e-12, 2000.0)
    ARR_ac(1.0e-12, -3.0)
    ARR_abc(1.0e-12, 2000.0, -3.0)
    """
    kpp_str = 'ARR_ab(1.0e-12, 2000.0)'
    logger.debug(kpp_str)
    kpp_to_micm.parse_kpp_arrhenius(kpp_str)
    assert True
