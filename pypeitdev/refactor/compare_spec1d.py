import numpy as np

from pypeit import specobjs
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]


sobjs1 = specobjs.SpecObjs.from_fitsfile(file1)
sobjs2 = specobjs.SpecObjs.from_fitsfile(file2)
if len(sys.argv) > 3:
    verbose=True
else:
    verbose=False

verbose_threshold = 1e-5

for sobj1 in sobjs1:
    idx = sobjs2.NAME == sobj1.NAME
    if np.sum(idx) == 0:
        print(f"{sobj1.NAME} is not in the second file")
        continue
    if np.sum(idx) > 1:
        print(f"{sobj1.NAME} is duplicated in the second file")
        continue

    sobj2 = sobjs2[idx][0]

    datamodel = sobj1.datamodel
    rtols = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000]
    for key in datamodel:
        if sobj1[key] is None or sobj2[key] is None:
            if sobj1[key] is None and sobj2[key] is None:
                if verbose:
                    print(f"{sobj1.NAME}.{key} are both None")
            elif sobj1[key] is None:
                print(f"{sobj1.NAME}.{key} differs, first file has None")
            else:
                print(f"{sobj1.NAME}.{key} differs, second file has None")
            continue

        if datamodel[key]["otype"]==np.ndarray:
            if datamodel[key]["atype"] == float:
                matched = False
                # Deal with None's by converting them to NaNs
                a1 = sobj1[key]
                a1[a1==None] = np.nan
                a2 = sobj2[key]
                a2[a2==None] = np.nan

                for rtol in rtols:
                    if np.allclose(a1, a2, rtol=rtol, equal_nan=True):
                        if verbose or rtol >= verbose_threshold:
                            print(f"{sobj1.NAME}.{key}: within {rtol}")
                        matched=True
                        break
                if not matched:
                    print(f"{sobj1.NAME}.{key}: differs by more than rtol=10000")
            else:
                if np.array_equal(sobj1[key], sobj2[key]):
                    if verbose:
                        print(f"{sobj1.NAME}.{key}: is equal")
                else:
                    print(f"{sobj1.NAME}.{key}: differs")
            
        elif datamodel[key]["otype"] == float:
            matched=False
            for rtol in rtols:
                if np.isclose(sobj1[key], sobj2[key], rtol=rtol):
                    if verbose or rtol >= verbose_threshold:
                        print(f"{sobj1.NAME}.{key}: within {rtol}")
                    matched=True
                    break
            if not matched:
                print(f"{sobj1.NAME}.{key}: differs by more than rtol=1")
        else:
            if sobj1[key] == sobj2[key]:
                if verbose:
                    print(f"{sobj1.NAME}.{key} is equal.")
            else:
                print(f"{sobj1.NAME}.{key} differs.")

            
for sobj2 in sobjs2:
    idx = sobjs1.NAME == sobj2.NAME 
    if np.sum(idx) == 0:
        print(f"{sobj2.NAME} is not in the first file")
    elif np.sum(idx) > 1:
        print(f"{sobj2.NAME} is duplicated in the first file")
