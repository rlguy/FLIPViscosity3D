#ifndef LEVELSETUTILS_H
#define LEVELSETUTILS_H

namespace LevelsetUtils {

extern float fractionInside(float phi_left, float phi_right);
extern float fractionInside(float phi_bl, float phi_br, float phi_tl, float phi_tr);

extern void _cycleArray(float* arr, int size);

}

#endif