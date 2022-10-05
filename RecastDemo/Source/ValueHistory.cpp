#include "ValueHistory.h"
#include "imgui.h"
#include <string.h>
#include <stdio.h>

#ifdef WIN32
#	define snprintf _snprintf
#endif

ValueHistory::ValueHistory() :
	m_hsamples(0)
{
	for (int i = 0; i < MAX_HISTORY; ++i)
		m_samples[i] = 0;
}

double ValueHistory::getSampleMin() const
{
	double val = m_samples[0];
	for (int i = 1; i < MAX_HISTORY; ++i)
		if (m_samples[i] < val)
			val = m_samples[i];
	return val;
} 

double ValueHistory::getSampleMax() const
{
	double val = m_samples[0];
	for (int i = 1; i < MAX_HISTORY; ++i)
		if (m_samples[i] > val)
			val = m_samples[i];
	return val;
}

double ValueHistory::getAverage() const
{
	double val = 0;
	for (int i = 0; i < MAX_HISTORY; ++i)
		val += m_samples[i];
	return val/(double)MAX_HISTORY;
}

void GraphParams::setRect(int ix, int iy, int iw, int ih, int ipad)
{
	x = ix;
	y = iy;
	w = iw;
	h = ih;
	pad = ipad;
}

void GraphParams::setValueRange(double ivmin, double ivmax, int indiv, const char* iunits)
{
	vmin = ivmin;
	vmax = ivmax;
	ndiv = indiv;
	strcpy(units, iunits);
}

void drawGraphBackground(const GraphParams* p)
{
	// BG
	imguiDrawRoundedRect((double)p->x, (double)p->y, (double)p->w, (double)p->h, (double)p->pad, imguiRGBA(64,64,64,128));
	
	const double sy = (p->h-p->pad*2) / (p->vmax-p->vmin);
	const double oy = p->y+p->pad-p->vmin*sy;
	
	char text[64];
	
	// Divider Lines
	for (int i = 0; i <= p->ndiv; ++i)
	{
		const double u = (double)i/(double)p->ndiv;
		const double v = p->vmin + (p->vmax-p->vmin)*u;
		snprintf(text, 64, "%.2lf %s", v, p->units);
		const double fy = oy + v*sy;
		imguiDrawText(p->x + p->w - p->pad, (int)fy-4, IMGUI_ALIGN_RIGHT, text, imguiRGBA(0,0,0,255));
		imguiDrawLine((double)p->x + (double)p->pad, fy, (double)p->x + (double)p->w - (double)p->pad - 50, fy, 1.0, imguiRGBA(0,0,0,64));
	}
}

void drawGraph(const GraphParams* p, const ValueHistory* graph,
			   int idx, const char* label, const unsigned int col)
{
	const double sx = (p->w - p->pad*2) / (double)graph->getSampleCount();
	const double sy = (p->h - p->pad*2) / (p->vmax - p->vmin);
	const double ox = (double)p->x + (double)p->pad;
	const double oy = (double)p->y + (double)p->pad - p->vmin*sy;
	
	// Values
	double px=0, py=0;
	for (int i = 0; i < graph->getSampleCount()-1; ++i)
	{
		const double x = ox + i*sx;
		const double y = oy + graph->getSample(i)*sy;
		if (i > 0)
			imguiDrawLine(px,py, x,y, 2.0, col);
		px = x;
		py = y;
	}
	
	// Label
	const int size = 15;
	const int spacing = 10;
	int ix = p->x + p->w + 5;
	int iy = p->y + p->h - (idx+1)*(size+spacing);
	
	imguiDrawRoundedRect((double)ix, (double)iy, (double)size, (double)size, 2.0, col);
	
	char text[64];
	snprintf(text, 64, "%.2f %s", graph->getAverage(), p->units);
	imguiDrawText(ix+size+5, iy+3, IMGUI_ALIGN_LEFT, label, imguiRGBA(255,255,255,192));
	imguiDrawText(ix+size+150, iy+3, IMGUI_ALIGN_RIGHT, text, imguiRGBA(255,255,255,128));
}

