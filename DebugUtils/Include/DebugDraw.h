//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef DEBUGDRAW_H
#define DEBUGDRAW_H

// Some math headers don't have PI defined.
static const double DU_PI = 3.141592653589793238462643383279502884;

enum duDebugDrawPrimitives
{
	DU_DRAW_POINTS,
	DU_DRAW_LINES,
	DU_DRAW_TRIS,
	DU_DRAW_QUADS,	
};

/// Abstract debug draw interface.
struct duDebugDraw
{
	virtual ~duDebugDraw() = 0;
	
	virtual void depthMask(bool state) = 0;

	virtual void texture(bool state) = 0;

	/// Begin drawing primitives.
	///  @param prim [in] primitive type to draw, one of rcDebugDrawPrimitives.
	///  @param size [in] size of a primitive, applies to point size and line width only.
	virtual void begin(duDebugDrawPrimitives prim, double size = 1.0) = 0;

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	virtual void vertex(const double* pos, unsigned int color) = 0;

	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	virtual void vertex(const double x, const double y, const double z, unsigned int color) = 0;

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	virtual void vertex(const double* pos, unsigned int color, const double* uv) = 0;
	
	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	virtual void vertex(const double x, const double y, const double z, unsigned int color, const double u, const double v) = 0;
	
	/// End drawing primitives.
	virtual void end() = 0;

	/// Compute a color for given area.
	virtual unsigned int areaToCol(unsigned int area);
};

inline unsigned int duRGBA(int r, int g, int b, int a)
{
	return ((unsigned int)r) | ((unsigned int)g << 8) | ((unsigned int)b << 16) | ((unsigned int)a << 24);
}

inline unsigned int duRGBAf(double fr, double fg, double fb, double fa)
{
	unsigned char r = (unsigned char)(fr*255.0);
	unsigned char g = (unsigned char)(fg*255.0);
	unsigned char b = (unsigned char)(fb*255.0);
	unsigned char a = (unsigned char)(fa*255.0);
	return duRGBA(r,g,b,a);
}

unsigned int duIntToCol(int i, int a);
void duIntToCol(int i, double* col);

inline unsigned int duMultCol(const unsigned int col, const unsigned int d)
{
	const unsigned int r = col & 0xff;
	const unsigned int g = (col >> 8) & 0xff;
	const unsigned int b = (col >> 16) & 0xff;
	const unsigned int a = (col >> 24) & 0xff;
	return duRGBA((r*d) >> 8, (g*d) >> 8, (b*d) >> 8, a);
}

inline unsigned int duDarkenCol(unsigned int col)
{
	return ((col >> 1) & 0x007f7f7f) | (col & 0xff000000);
}

inline unsigned int duLerpCol(unsigned int ca, unsigned int cb, unsigned int u)
{
	const unsigned int ra = ca & 0xff;
	const unsigned int ga = (ca >> 8) & 0xff;
	const unsigned int ba = (ca >> 16) & 0xff;
	const unsigned int aa = (ca >> 24) & 0xff;
	const unsigned int rb = cb & 0xff;
	const unsigned int gb = (cb >> 8) & 0xff;
	const unsigned int bb = (cb >> 16) & 0xff;
	const unsigned int ab = (cb >> 24) & 0xff;
	
	unsigned int r = (ra*(255-u) + rb*u)/255;
	unsigned int g = (ga*(255-u) + gb*u)/255;
	unsigned int b = (ba*(255-u) + bb*u)/255;
	unsigned int a = (aa*(255-u) + ab*u)/255;
	return duRGBA(r,g,b,a);
}

inline unsigned int duTransCol(unsigned int c, unsigned int a)
{
	return (a<<24) | (c & 0x00ffffff);
}


void duCalcBoxColors(unsigned int* colors, unsigned int colTop, unsigned int colSide);

void duDebugDrawCylinderWire(struct duDebugDraw* dd, double minx, double miny, double minz,
							 double maxx, double maxy, double maxz, unsigned int col, const double lineWidth);

void duDebugDrawBoxWire(struct duDebugDraw* dd, double minx, double miny, double minz,
						double maxx, double maxy, double maxz, unsigned int col, const double lineWidth);

void duDebugDrawArc(struct duDebugDraw* dd, const double x0, const double y0, const double z0,
					const double x1, const double y1, const double z1, const double h,
					const double as0, const double as1, unsigned int col, const double lineWidth);

void duDebugDrawArrow(struct duDebugDraw* dd, const double x0, const double y0, const double z0,
					  const double x1, const double y1, const double z1,
					  const double as0, const double as1, unsigned int col, const double lineWidth);

void duDebugDrawCircle(struct duDebugDraw* dd, const double x, const double y, const double z,
					   const double r, unsigned int col, const double lineWidth);

void duDebugDrawCross(struct duDebugDraw* dd, const double x, const double y, const double z,
					  const double size, unsigned int col, const double lineWidth);

void duDebugDrawBox(struct duDebugDraw* dd, double minx, double miny, double minz,
					double maxx, double maxy, double maxz, const unsigned int* fcol);

void duDebugDrawCylinder(struct duDebugDraw* dd, double minx, double miny, double minz,
						 double maxx, double maxy, double maxz, unsigned int col);

void duDebugDrawGridXZ(struct duDebugDraw* dd, const double ox, const double oy, const double oz,
					   const int w, const int h, const double size,
					   const unsigned int col, const double lineWidth);


// Versions without begin/end, can be used to draw multiple primitives.
void duAppendCylinderWire(struct duDebugDraw* dd, double minx, double miny, double minz,
						  double maxx, double maxy, double maxz, unsigned int col);

void duAppendBoxWire(struct duDebugDraw* dd, double minx, double miny, double minz,
					 double maxx, double maxy, double maxz, unsigned int col);

void duAppendBoxPoints(struct duDebugDraw* dd, double minx, double miny, double minz,
					   double maxx, double maxy, double maxz, unsigned int col);

void duAppendArc(struct duDebugDraw* dd, const double x0, const double y0, const double z0,
				 const double x1, const double y1, const double z1, const double h,
				 const double as0, const double as1, unsigned int col);

void duAppendArrow(struct duDebugDraw* dd, const double x0, const double y0, const double z0,
				   const double x1, const double y1, const double z1,
				   const double as0, const double as1, unsigned int col);

void duAppendCircle(struct duDebugDraw* dd, const double x, const double y, const double z,
					const double r, unsigned int col);

void duAppendCross(struct duDebugDraw* dd, const double x, const double y, const double z,
				   const double size, unsigned int col);

void duAppendBox(struct duDebugDraw* dd, double minx, double miny, double minz,
				 double maxx, double maxy, double maxz, const unsigned int* fcol);

void duAppendCylinder(struct duDebugDraw* dd, double minx, double miny, double minz,
					  double maxx, double maxy, double maxz, unsigned int col);


class duDisplayList : public duDebugDraw
{
	double* m_pos;
	unsigned int* m_color;
	int m_size;
	int m_cap;

	bool m_depthMask;
	duDebugDrawPrimitives m_prim;
	double m_primSize;
	
	void resize(int cap);
	
public:
	duDisplayList(int cap = 512);
	~duDisplayList();
	virtual void depthMask(bool state);
	virtual void begin(duDebugDrawPrimitives prim, double size = 1.0);
	virtual void vertex(const double x, const double y, const double z, unsigned int color);
	virtual void vertex(const double* pos, unsigned int color);
	virtual void end();
	void clear();
	void draw(struct duDebugDraw* dd);
private:
	// Explicitly disabled copy constructor and copy assignment operator.
	duDisplayList(const duDisplayList&);
	duDisplayList& operator=(const duDisplayList&);
};


#endif // DEBUGDRAW_H
