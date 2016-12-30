/*
    LevelSet methods addapted from Christopher Batty's levelset_util.cpp:
        https://github.com/christopherbatty/Fluid3D/blob/master/levelset_util.cpp

*/

#include "levelsetutils.h"

namespace LevelsetUtils {

//Given two signed distance values (line endpoints), determine what fraction of a connecting segment is "inside"
float fractionInside(float phiLeft, float phiRight) {
    if(phiLeft < 0 && phiRight < 0) {
        return 1;
    }
    if (phiLeft < 0 && phiRight >= 0) {
        return phiLeft / (phiLeft - phiRight);
    }
    if(phiLeft >= 0 && phiRight < 0) {
        return phiRight / (phiRight - phiLeft);
    }
        
    return 0;
}

void _cycleArray(float* arr, int size) {
    float t = arr[0];
    for(int i = 0; i < size - 1; ++i) {
        arr[i] = arr[i + 1];
    }
    arr[size - 1] = t;
}

//Given four signed distance values (square corners), determine what fraction of the square is "inside"
float fractionInside(float phibl, float phibr, float phitl, float phitr) {
    
    int insideCount = (phibl < 0 ? 1 : 0) + 
                      (phitl < 0 ? 1 : 0) + 
                      (phibr < 0 ? 1 : 0) + 
                      (phitr < 0 ? 1 : 0);
    float list[] = { phibl, phibr, phitr, phitl };

    if(insideCount == 4) {
        return 1;
    } else if (insideCount == 3) {
        //rotate until the positive value is in the first position
        while(list[0] < 0) {
            _cycleArray(list, 4);
        }

        //Work out the area of the exterior triangle
        float side0 = 1 - fractionInside(list[0], list[3]);
        float side1 = 1 - fractionInside(list[0], list[1]);
        return 1.0f - 0.5f * side0 * side1;
    } else if(insideCount == 2) {
        
        //rotate until a negative value is in the first position, and the next negative is in either slot 1 or 2.
        while(list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
            _cycleArray(list , 4);
        } 
        
        if(list[1] < 0) { //the matching signs are adjacent
            float sideLeft = fractionInside(list[0], list[3]);
            float sideRight = fractionInside(list[1], list[2]);
            return  0.5f * (sideLeft + sideRight);
        } else { 
            //matching signs are diagonally opposite
            //determine the centre point's sign to disambiguate this case
            float middlePoint = 0.25f * (list[0] + list[1] + list[2] + list[3]);
            if(middlePoint < 0) {
                float area = 0;

                //first triangle (top left)
                float side1 = 1 - fractionInside(list[0], list[3]);
                float side3 = 1 - fractionInside(list[2], list[3]);

                area += 0.5f * side1 * side3;

                //second triangle (top right)
                float side2 = 1 - fractionInside(list[2], list[1]);
                float side0 = 1 - fractionInside(list[0], list[1]);
                area += 0.5f * side0 * side2;
                
                return 1.0f - area;
            }
            else {
                float area = 0;

                //first triangle (bottom left)
                float side0 = fractionInside(list[0], list[1]);
                float side1 = fractionInside(list[0], list[3]);
                area += 0.5f * side0*side1;

                //second triangle (top right)
                float side2 = fractionInside(list[2], list[1]);
                float side3 = fractionInside(list[2], list[3]);
                area += 0.5f * side2 * side3;
                return area;
            }
            
        }
    } else if(insideCount == 1) {
        //rotate until the negative value is in the first position
        while(list[0] >= 0) {
            _cycleArray(list, 4);
        }

        //Work out the area of the interior triangle, and subtract from 1.
        float side0 = fractionInside(list[0], list[3]);
        float side1 = fractionInside(list[0], list[1]);
        return 0.5f * side0 * side1;
    } else {
        return 0;
    }

}

}