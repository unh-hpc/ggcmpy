"""
Download data from CDAWeb
"""
from __future__ import absolute_import, print_function
import os as _os

from . import low_level as _ll


def _getData(starttime, endtime, directory, mapping):
    url = 'http://cdaweb.gsfc.nasa.gov/WS/cdasr/1'
    dv = 'sp_phys'
    files = []
    for dataset, prefix in mapping.items():
        files.extend(_ll.downloadResults(_ll.getTextData(url, dv, starttime,
                                                         endtime, dataset),
                                         prefix=prefix, directory=directory))
    return files

def getAceData(starttime, endtime, directory=_os.curdir):
    """
    This routine should get the ACE data suitable
    for an OpenGGCM run
    """
    d = {"AC_H0_MFI/BGSEc":"B_",
         "AC_H0_MFI/SC_pos_GSE":"Position_",
         "AC_H0_SWE/V_GSE":"V_",
         "AC_H0_SWE/Tpr":"T_",
         "AC_H0_SWE/Np":"N_",
         }
    return _getData(starttime, endtime, directory, d)

def getWindData(starttime, endtime, directory=_os.curdir):
    """
    This routine should get the WIND data suitable
    for an OpenGGCM run
    """
    d = {"WI_H0_MFI/B3GSE":"B_",
         "WI_H0_MFI/P1GSE":"Position_",
         "WI_K0_SWE/V_GSE":"V_",
         "WI_K0_SWE/THERMAL_SPD":"Vth_",
         "WI_K0_SWE/Np":"N_",
         }
    return _getData(starttime, endtime, directory, d)


def getGeotailData(starttime, endtime, directory=_os.curdir):
    d = {'GE_H0_CPI/SW_Vc':"V_",
         'GE_H0_CPI/SW_P_Den':"N_",
         'GE_H0_CPI/SW_T':'T_',
         'GE_H0_CPI/SW_Pressure':'pp_',
         'GE_K0_MGF/IB_vector':'B_',
         'GE_K0_MGF/POS':"Position_"}
    return _getData(starttime, endtime, directory, d)


def getOMNIData(starttime, endtime, directory=_os.curdir):
    """
    This routine should get the OMNI data suitable
    for an OpenGGCM run
    """
    d = {"OMNI_HRO_1MIN/T":"T_",
         "OMNI_HRO_1MIN/Pressure":"pp_",
         "OMNI_HRO_1MIN/proton_density":"rr_",
         "OMNI_HRO_1MIN/Vx":"Vx_",
         "OMNI_HRO_1MIN/Vy":"Vy_",
         "OMNI_HRO_1MIN/Vz":"Vz_",
         "OMNI_HRO_1MIN/BX_GSE":"Bx_",
         "OMNI_HRO_1MIN/BY_GSE":"By_",
         "OMNI_HRO_1MIN/BZ_GSE":"Bz_",
          }
    return _getData(starttime, endtime, directory, d)

def getAuroralIndicesData(starttime, endtime, directory=_os.curdir):
    """
    This routine should get the AE, AL,  and AU indices
    from WDC Kyoto (via OMNI)
    """
    d = {"OMNI_HRO_1MIN/AU_INDEX":"AU_",
         "OMNI_HRO_1MIN/AL_INDEX":"AL_",
         "OMNI_HRO_1MIN/AE_INDEX":"AE_",
         }
    return _getData(starttime, endtime, directory, d)

def getSymIndicesData(starttime, endtime, directory=_os.curdir):
    """
    This routine should get the SYM/ASY-H/D indices
    from WDC Kyoto (via OMNI)
    """
    d = {"OMNI_HRO_1MIN/SYM_D":"SYMD_",
         "OMNI_HRO_1MIN/SYM_H":"SYMH_",
         "OMNI_HRO_1MIN/ASY_D":"ASYD_",
         "OMNI_HRO_1MIN/ASY_H":"ASYH_",
         "OMNI2_H0_MRG1HR/DST1800":"DST_",
         }
    return _getData(starttime, endtime, directory, d)
