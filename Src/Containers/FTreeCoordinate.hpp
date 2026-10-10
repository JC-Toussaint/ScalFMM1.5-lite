// See LICENCE file at project root
#ifndef FTREECOORDINATE_HPP
#define FTREECOORDINATE_HPP

#include <array>

#include "../Utils/FGlobal.hpp"

/**
 * @class FTreeCoordinate
 * Coordinates (x, y, z) of a cell in the octree grid of its level, and the conversions
 * to and from its Morton index (bits of x, y, z interleaved, x being the most significant).
 */
class FTreeCoordinate {
    std::array<int, 3> coord{};   //< x, y, z

public:
    FTreeCoordinate() = default;

    /** From a Morton index */
    explicit FTreeCoordinate(const MortonIndex mindex) {
        setPositionFromMorton(mindex);
    }

    explicit FTreeCoordinate(const int inX, const int inY, const int inZ) : coord{inX, inY, inZ} {}

    int getX() const { return coord[0]; }
    int getY() const { return coord[1]; }
    int getZ() const { return coord[2]; }

    void setX(const int inX) { coord[0] = inX; }
    void setY(const int inY) { coord[1] = inY; }
    void setZ(const int inZ) { coord[2] = inZ; }

    /** Morton index of the coordinates */
    MortonIndex getMortonIndex() const {
        MortonIndex index = 0x0LL;
        MortonIndex mask = 0x1LL;
        MortonIndex mx = MortonIndex(coord[0]) << 2;
        MortonIndex my = MortonIndex(coord[1]) << 1;
        MortonIndex mz = coord[2];
        while( (mask <= mz)
               || ((mask << 1) <= my)
               || ((mask << 2) <= mx))
        {
            index |= (mz & mask);
            mask <<= 1;
            index |= (my & mask);
            mask <<= 1;
            index |= (mx & mask);
            mask <<= 1;
            mz <<= 2;
            my <<= 2;
            mx <<= 2;
        }
        return index;
    }

    /** Coordinates from a Morton index */
    void setPositionFromMorton(MortonIndex inIndex) {
        MortonIndex mask = 0x1LL;
        coord = {0, 0, 0};
        while(inIndex >= mask) {
            coord[2] |= int(inIndex & mask);
            inIndex >>= 1;
            coord[1] |= int(inIndex & mask);
            inIndex >>= 1;
            coord[0] |= int(inIndex & mask);
            mask <<= 1;
        }
    }
};

#endif // FTREECOORDINATE_HPP
