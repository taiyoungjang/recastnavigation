#include "catch.hpp"

#include "DetourCommon.h"

TEST_CASE("dtRandomPointInConvexPoly")
{
	SECTION("Properly works when the argument 's' is 1.0")
	{
		const double pts[] = {
			0, 0, 0,
			0, 0, 1,
			1, 0, 0,
		};
		const int npts = 3;
		double areas[6];
		double out[3];

		dtRandomPointInConvexPoly(pts, npts, areas, 0.0, 1.0, out);
		REQUIRE(out[0] == Approx(0));
		REQUIRE(out[1] == Approx(0));
		REQUIRE(out[2] == Approx(1));

		dtRandomPointInConvexPoly(pts, npts, areas, 0.5, 1.0, out);
		REQUIRE(out[0] == Approx(1.0 / 2));
		REQUIRE(out[1] == Approx(0));
		REQUIRE(out[2] == Approx(1.0 / 2));

		dtRandomPointInConvexPoly(pts, npts, areas, 1.0, 1.0, out);
		REQUIRE(out[0] == Approx(1));
		REQUIRE(out[1] == Approx(0));
		REQUIRE(out[2] == Approx(0));
	}
}
