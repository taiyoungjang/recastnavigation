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

#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>

#include "SDL.h"
#include "SDL_opengl.h"
#ifdef __APPLE__
#	include <OpenGL/glu.h>
#else
#	include <GL/glu.h>
#endif

#include <vector>
#include <string>

#include "imgui.h"
#include "imguiRenderGL.h"

#include "Recast.h"
#include "RecastDebugDraw.h"
#include "InputGeom.h"
#include "TestCase.h"
#include "Filelist.h"
#include "Sample_SoloMesh.h"
#include "Sample_TileMesh.h"
#include "Sample_TempObstacles.h"
#include "Sample_Debug.h"

#ifdef WIN32
#	define snprintf _snprintf
#	define putenv _putenv
#endif

using std::string;
using std::vector;

struct SampleItem
{
	Sample* (*create)();
	const string name;
};
Sample* createSolo() { return new Sample_SoloMesh(); }
Sample* createTile() { return new Sample_TileMesh(); }
Sample* createTempObstacle() { return new Sample_TempObstacles(); }
Sample* createDebug() { return new Sample_Debug(); }
static SampleItem g_samples[] =
{
	{ createSolo, "Solo Mesh" },
	{ createTile, "Tile Mesh" },
	{ createTempObstacle, "Temp Obstacles" },
};
static const int g_nsamples = sizeof(g_samples) / sizeof(SampleItem);

int main(int /*argc*/, char** /*argv*/)
{
	// Init SDL
	if (SDL_Init(SDL_INIT_EVERYTHING) != 0)
	{
		printf("Could not initialise SDL.\nError: %s\n", SDL_GetError());
		return -1;
	}

    // Use OpenGL render driver.
    SDL_SetHint(SDL_HINT_RENDER_DRIVER, "opengl");

	// Enable depth buffer.
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	
	// Set color channel depth.
	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
	
	// 4x MSAA.
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);

	SDL_DisplayMode displayMode;
	SDL_GetCurrentDisplayMode(0, &displayMode);

	bool presentationMode = false;
	Uint32 flags = SDL_WINDOW_OPENGL;
	int width;
	int height;
	if (presentationMode)
	{
		// Create a fullscreen window at the native resolution.
		width = displayMode.w;
		height = displayMode.h;
		flags |= SDL_WINDOW_FULLSCREEN;
	}
	else
	{
		double aspect = 16.0 / 9.0;
		width = rcMin(displayMode.w, (int)(displayMode.h * aspect)) - 80;
		height = displayMode.h - 80;
	}
	
	SDL_Window* window;
	SDL_Renderer* renderer;
	int errorCode = SDL_CreateWindowAndRenderer(width, height, flags, &window, &renderer);

	if (errorCode != 0 || !window || !renderer)
	{
		printf("Could not initialise SDL opengl\nError: %s\n", SDL_GetError());
		return -1;
	}

	SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
	SDL_GL_CreateContext(window);

	if (!imguiRenderGLInit("DroidSans.ttf"))
	{
		printf("Could not init GUI renderer.\n");
		SDL_Quit();
		return -1;
	}
	
	double t = 0.0;
	double timeAcc = 0.0;
	Uint32 prevFrameTime = SDL_GetTicks();
	int mousePos[2] = {0, 0};
	int origMousePos[2] = {0, 0}; // Used to compute mouse movement totals across frames.
	
	double cameraEulers[] = {45, -45};
	double cameraPos[] = {0, 0, 0};
	double camr = 1000;
	double origCameraEulers[] = {0, 0}; // Used to compute rotational changes across frames.
	
	double moveFront = 0.0, moveBack = 0.0, moveLeft = 0.0, moveRight = 0.0, moveUp = 0.0, moveDown = 0.0;
	
	double scrollZoom = 0;
	bool rotate = false;
	bool movedDuringRotate = false;
	double rayStart[3];
	double rayEnd[3];
	bool mouseOverMenu = false;
	
	bool showMenu = !presentationMode;
	bool showLog = false;
	bool showTools = true;
	bool showLevels = false;
	bool showSample = false;
	bool showTestCases = false;

	// Window scroll positions.
	int propScroll = 0;
	int logScroll = 0;
	int toolsScroll = 0;
	
	string sampleName = "Choose Sample...";
	
	vector<string> files;
	const string meshesFolder = "Meshes";
	string meshName = "Choose Mesh...";
	
	double markerPosition[3] = {0, 0, 0};
	bool markerPositionSet = false;
	
	InputGeom* geom = 0;
	Sample* sample = 0;

	const string testCasesFolder = "TestCases";
	TestCase* test = 0;

	BuildContext ctx;
	
	// Fog.
	float fogColor[4] = { 0.32, 0.31, 0.30, 1.0 };
	glEnable(GL_FOG);
	glFogi(GL_FOG_MODE, GL_LINEAR);
	glFogf(GL_FOG_START, camr * 0.1);
	glFogf(GL_FOG_END, camr * 1.25);
	glFogfv(GL_FOG_COLOR, fogColor);
	
	glEnable(GL_CULL_FACE);
	glDepthFunc(GL_LEQUAL);
	
	bool done = false;
	while(!done)
	{
		// Handle input events.
		int mouseScroll = 0;
		bool processHitTest = false;
		bool processHitTestShift = false;
		SDL_Event event;
		
		while (SDL_PollEvent(&event))
		{
			switch (event.type)
			{
				case SDL_KEYDOWN:
					// Handle any key presses here.
					if (event.key.keysym.sym == SDLK_ESCAPE)
					{
						done = true;
					}
					else if (event.key.keysym.sym == SDLK_t)
					{
						showLevels = false;
						showSample = false;
						showTestCases = true;
						scanDirectory(testCasesFolder, ".txt", files);
					}
					else if (event.key.keysym.sym == SDLK_TAB)
					{
						showMenu = !showMenu;
					}
					else if (event.key.keysym.sym == SDLK_SPACE)
					{
						if (sample)
							sample->handleToggle();
					}
					else if (event.key.keysym.sym == SDLK_1)
					{
						if (sample)
							sample->handleStep();
					}
					else if (event.key.keysym.sym == SDLK_9)
					{
						if (sample && geom)
						{
							string savePath = meshesFolder + "/";
							BuildSettings settings;
							memset(&settings, 0, sizeof(settings));

							rcVcopy(settings.navMeshBMin, geom->getNavMeshBoundsMin());
							rcVcopy(settings.navMeshBMax, geom->getNavMeshBoundsMax());

							sample->collectSettings(settings);

							geom->saveGeomSet(&settings);
						}
					}
					break;
				
				case SDL_MOUSEWHEEL:
					if (event.wheel.y < 0)
					{
						// wheel down
						if (mouseOverMenu)
						{
							mouseScroll++;
						}
						else
						{
							scrollZoom += 1.0;
						}
					}
					else
					{
						if (mouseOverMenu)
						{
							mouseScroll--;
						}
						else
						{
							scrollZoom -= 1.0;
						}
					}
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_RIGHT)
					{
						if (!mouseOverMenu)
						{
							// Rotate view
							rotate = true;
							movedDuringRotate = false;
							origMousePos[0] = mousePos[0];
							origMousePos[1] = mousePos[1];
							origCameraEulers[0] = cameraEulers[0];
							origCameraEulers[1] = cameraEulers[1];
						}
					}
					break;
					
				case SDL_MOUSEBUTTONUP:
					// Handle mouse clicks here.
					if (event.button.button == SDL_BUTTON_RIGHT)
					{
						rotate = false;
						if (!mouseOverMenu)
						{
							if (!movedDuringRotate)
							{
								processHitTest = true;
								processHitTestShift = true;
							}
						}
					}
					else if (event.button.button == SDL_BUTTON_LEFT)
					{
						if (!mouseOverMenu)
						{
							processHitTest = true;
							processHitTestShift = (SDL_GetModState() & KMOD_SHIFT) ? true : false;
						}
					}
					
					break;
					
				case SDL_MOUSEMOTION:
					mousePos[0] = event.motion.x;
					mousePos[1] = height-1 - event.motion.y;
					
					if (rotate)
					{
						int dx = mousePos[0] - origMousePos[0];
						int dy = mousePos[1] - origMousePos[1];
						cameraEulers[0] = origCameraEulers[0] - dy * 0.25;
						cameraEulers[1] = origCameraEulers[1] + dx * 0.25;
						if (dx * dx + dy * dy > 3 * 3)
						{
							movedDuringRotate = true;
						}
					}
					break;
					
				case SDL_QUIT:
					done = true;
					break;
					
				default:
					break;
			}
		}

		unsigned char mouseButtonMask = 0;
		if (SDL_GetMouseState(0, 0) & SDL_BUTTON_LMASK)
			mouseButtonMask |= IMGUI_MBUT_LEFT;
		if (SDL_GetMouseState(0, 0) & SDL_BUTTON_RMASK)
			mouseButtonMask |= IMGUI_MBUT_RIGHT;
		
		Uint32 time = SDL_GetTicks();
		double dt = (time - prevFrameTime) / 1000.0;
		prevFrameTime = time;
		
		t += dt;

		// Hit test mesh.
		if (processHitTest && geom && sample)
		{
			double hitTime;
			bool hit = geom->raycastMesh(rayStart, rayEnd, hitTime);
			
			if (hit)
			{
				if (SDL_GetModState() & KMOD_CTRL)
				{
					// Marker
					markerPositionSet = true;
					markerPosition[0] = rayStart[0] + (rayEnd[0] - rayStart[0]) * hitTime;
					markerPosition[1] = rayStart[1] + (rayEnd[1] - rayStart[1]) * hitTime;
					markerPosition[2] = rayStart[2] + (rayEnd[2] - rayStart[2]) * hitTime;
				}
				else
				{
					double pos[3];
					pos[0] = rayStart[0] + (rayEnd[0] - rayStart[0]) * hitTime;
					pos[1] = rayStart[1] + (rayEnd[1] - rayStart[1]) * hitTime;
					pos[2] = rayStart[2] + (rayEnd[2] - rayStart[2]) * hitTime;
					sample->handleClick(rayStart, pos, processHitTestShift);
				}
			}
			else
			{
				if (SDL_GetModState() & KMOD_CTRL)
				{
					// Marker
					markerPositionSet = false;
				}
			}
		}
		
		// Update sample simulation.
		const double SIM_RATE = 20;
		const double DELTA_TIME = 1.0 / SIM_RATE;
        const double one = 1.0;
		timeAcc = rcClamp(timeAcc + dt, -one, one);
		int simIter = 0;
		while (timeAcc > DELTA_TIME)
		{
			timeAcc -= DELTA_TIME;
			if (simIter < 5 && sample)
			{
				sample->handleUpdate(DELTA_TIME);
			}
			simIter++;
		}

		// Clamp the framerate so that we do not hog all the CPU.
		const double MIN_FRAME_TIME = 1.0 / 40.0;
		if (dt < MIN_FRAME_TIME)
		{
			int ms = (int)((MIN_FRAME_TIME - dt) * 1000.0);
			if (ms > 10) ms = 10;
			if (ms >= 0) SDL_Delay(ms);
		}
		
		// Set the viewport.
		glViewport(0, 0, width, height);
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		
		// Clear the screen
		glClearColor(0.3, 0.3, 0.32, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_DEPTH_TEST);
		
		// Compute the projection matrix.
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(50.0, (double)width/(double)height, 1.0, camr);
		GLdouble projectionMatrix[16];
		glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);
		
		// Compute the modelview matrix.
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(cameraEulers[0], 1, 0, 0);
		glRotatef(cameraEulers[1], 0, 1, 0);
		glTranslatef(-cameraPos[0], -cameraPos[1], -cameraPos[2]);
		GLdouble modelviewMatrix[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMatrix);
		
		// Get hit ray position and direction.
		GLdouble x, y, z;
		gluUnProject(mousePos[0], mousePos[1], 0.0, modelviewMatrix, projectionMatrix, viewport, &x, &y, &z);
		rayStart[0] = (double)x;
		rayStart[1] = (double)y;
		rayStart[2] = (double)z;
		gluUnProject(mousePos[0], mousePos[1], 1.0, modelviewMatrix, projectionMatrix, viewport, &x, &y, &z);
		rayEnd[0] = (double)x;
		rayEnd[1] = (double)y;
		rayEnd[2] = (double)z;
		
		// Handle keyboard movement.
		const Uint8* keystate = SDL_GetKeyboardState(NULL);
        const double zero = 0.0;
		moveFront	= rcClamp(moveFront	+ dt * 4 * ((keystate[SDL_SCANCODE_W] || keystate[SDL_SCANCODE_UP		]) ? 1 : -1), zero, one);
		moveLeft	= rcClamp(moveLeft	+ dt * 4 * ((keystate[SDL_SCANCODE_A] || keystate[SDL_SCANCODE_LEFT		]) ? 1 : -1), zero, one);
		moveBack	= rcClamp(moveBack	+ dt * 4 * ((keystate[SDL_SCANCODE_S] || keystate[SDL_SCANCODE_DOWN		]) ? 1 : -1), zero, one);
		moveRight	= rcClamp(moveRight	+ dt * 4 * ((keystate[SDL_SCANCODE_D] || keystate[SDL_SCANCODE_RIGHT	]) ? 1 : -1), zero, one);
		moveUp		= rcClamp(moveUp	+ dt * 4 * ((keystate[SDL_SCANCODE_Q] || keystate[SDL_SCANCODE_PAGEUP	]) ? 1 : -1), zero, one);
		moveDown	= rcClamp(moveDown	+ dt * 4 * ((keystate[SDL_SCANCODE_E] || keystate[SDL_SCANCODE_PAGEDOWN	]) ? 1 : -1), zero, one);
		
		double keybSpeed = 22.0;
		if (SDL_GetModState() & KMOD_SHIFT)
		{
			keybSpeed *= 4.0;
		}
		
		double movex = (moveRight - moveLeft) * keybSpeed * dt;
		double movey = (moveBack - moveFront) * keybSpeed * dt + scrollZoom * 2.0;
		scrollZoom = 0;
		
		cameraPos[0] += movex * (double)modelviewMatrix[0];
		cameraPos[1] += movex * (double)modelviewMatrix[4];
		cameraPos[2] += movex * (double)modelviewMatrix[8];
		
		cameraPos[0] += movey * (double)modelviewMatrix[2];
		cameraPos[1] += movey * (double)modelviewMatrix[6];
		cameraPos[2] += movey * (double)modelviewMatrix[10];

		cameraPos[1] += (moveUp - moveDown) * keybSpeed * dt;

		glEnable(GL_FOG);

		if (sample)
			sample->handleRender();
		if (test)
			test->handleRender();
		
		glDisable(GL_FOG);
		
		// Render GUI
		glDisable(GL_DEPTH_TEST);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(0, width, 0, height);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		mouseOverMenu = false;
		
		imguiBeginFrame(mousePos[0], mousePos[1], mouseButtonMask, mouseScroll);
		
		if (sample)
		{
			sample->handleRenderOverlay((double*)projectionMatrix, (double*)modelviewMatrix, (int*)viewport);
		}
		if (test)
		{
			if (test->handleRenderOverlay((double*)projectionMatrix, (double*)modelviewMatrix, (int*)viewport))
				mouseOverMenu = true;
		}

		// Help text.
		if (showMenu)
		{
			const char msg[] = "W/S/A/D: Move  RMB: Rotate";
			imguiDrawText(280, height-20, IMGUI_ALIGN_LEFT, msg, imguiRGBA(255,255,255,128));
		}
		
		if (showMenu)
		{
			if (imguiBeginScrollArea("Properties", width-250-10, 10, 250, height-20, &propScroll))
				mouseOverMenu = true;

			if (imguiCheck("Show Log", showLog))
				showLog = !showLog;
			if (imguiCheck("Show Tools", showTools))
				showTools = !showTools;

			imguiSeparator();
			imguiLabel("Sample");
			if (imguiButton(sampleName.c_str()))
			{
				if (showSample)
				{
					showSample = false;
				}
				else
				{
					showSample = true;
					showLevels = false;
					showTestCases = false;
				}
			}
			
			imguiSeparator();
			imguiLabel("Input Mesh");
			if (imguiButton(meshName.c_str()))
			{
				if (showLevels)
				{
					showLevels = false;
				}
				else
				{
					showSample = false;
					showTestCases = false;
					showLevels = true;
					scanDirectory(meshesFolder, ".obj", files);
					scanDirectoryAppend(meshesFolder, ".gset", files);
				}
			}
			if (geom)
			{
				char text[64];
				snprintf(text, 64, "Verts: %.1lfk  Tris: %.1lfk",
						 geom->getMesh()->getVertCount()/1000.0,
						 geom->getMesh()->getTriCount()/1000.0);
				imguiValue(text);
			}
			imguiSeparator();

			if (geom && sample)
			{
				imguiSeparatorLine();
				
				sample->handleSettings();

				if (imguiButton("Build"))
				{
					ctx.resetLog();
					if (!sample->handleBuild())
					{
						showLog = true;
						logScroll = 0;
					}
					ctx.dumpLog("Build log %s:", meshName.c_str());
					
					// Clear test.
					delete test;
					test = 0;
				}

				imguiSeparator();
			}
			
			if (sample)
			{
				imguiSeparatorLine();
				sample->handleDebugMode();
			}

			imguiEndScrollArea();
		}
		
		// Sample selection dialog.
		if (showSample)
		{
			static int levelScroll = 0;
			if (imguiBeginScrollArea("Choose Sample", width-10-250-10-200, height-10-250, 200, 250, &levelScroll))
				mouseOverMenu = true;

			Sample* newSample = 0;
			for (int i = 0; i < g_nsamples; ++i)
			{
				if (imguiItem(g_samples[i].name.c_str()))
				{
					newSample = g_samples[i].create();
					if (newSample)
						sampleName = g_samples[i].name;
				}
			}
			if (newSample)
			{
				delete sample;
				sample = newSample;
				sample->setContext(&ctx);
				if (geom)
				{
					sample->handleMeshChanged(geom);
				}
				showSample = false;
			}

			if (geom || sample)
			{
				const double* bmin = 0;
				const double* bmax = 0;
				if (geom)
				{
					bmin = geom->getNavMeshBoundsMin();
					bmax = geom->getNavMeshBoundsMax();
				}
				// Reset camera and fog to match the mesh bounds.
				if (bmin && bmax)
				{
					camr = sqrt(rcSqr(bmax[0]-bmin[0]) +
								 rcSqr(bmax[1]-bmin[1]) +
								 rcSqr(bmax[2]-bmin[2])) / 2;
					cameraPos[0] = (bmax[0] + bmin[0]) / 2 + camr;
					cameraPos[1] = (bmax[1] + bmin[1]) / 2 + camr;
					cameraPos[2] = (bmax[2] + bmin[2]) / 2 + camr;
					camr *= 3;
				}
				cameraEulers[0] = 45;
				cameraEulers[1] = -45;
				glFogf(GL_FOG_START, camr*0.1);
				glFogf(GL_FOG_END, camr*1.25);
			}
			
			imguiEndScrollArea();
		}
		
		// Level selection dialog.
		if (showLevels)
		{
			static int levelScroll = 0;
			if (imguiBeginScrollArea("Choose Level", width - 10 - 250 - 10 - 200, height - 10 - 450, 200, 450, &levelScroll))
				mouseOverMenu = true;
			
			vector<string>::const_iterator fileIter = files.begin();
			vector<string>::const_iterator filesEnd = files.end();
			vector<string>::const_iterator levelToLoad = filesEnd;
			for (; fileIter != filesEnd; ++fileIter)
			{
				if (imguiItem(fileIter->c_str()))
				{
					levelToLoad = fileIter;
				}
			}
			
			if (levelToLoad != filesEnd)
			{
				meshName = *levelToLoad;
				showLevels = false;
				
				delete geom;
				geom = 0;
				
				string path = meshesFolder + "/" + meshName;
				
				geom = new InputGeom;
				if (!geom->load(&ctx, path))
				{
					delete geom;
					geom = 0;

					// Destroy the sample if it already had geometry loaded, as we've just deleted it!
					if (sample && sample->getInputGeom())
					{
						delete sample;
						sample = 0;
					}
					
					showLog = true;
					logScroll = 0;
					ctx.dumpLog("Geom load log %s:", meshName.c_str());
				}
				if (sample && geom)
				{
					sample->handleMeshChanged(geom);
				}

				if (geom || sample)
				{
					const double* bmin = 0;
					const double* bmax = 0;
					if (geom)
					{
						bmin = geom->getNavMeshBoundsMin();
						bmax = geom->getNavMeshBoundsMax();
					}
					// Reset camera and fog to match the mesh bounds.
					if (bmin && bmax)
					{
						camr = sqrt(rcSqr(bmax[0]-bmin[0]) +
									 rcSqr(bmax[1]-bmin[1]) +
									 rcSqr(bmax[2]-bmin[2])) / 2;
						cameraPos[0] = (bmax[0] + bmin[0]) / 2 + camr;
						cameraPos[1] = (bmax[1] + bmin[1]) / 2 + camr;
						cameraPos[2] = (bmax[2] + bmin[2]) / 2 + camr;
						camr *= 3;
					}
					cameraEulers[0] = 45;
					cameraEulers[1] = -45;
					glFogf(GL_FOG_START, camr * 0.1);
					glFogf(GL_FOG_END, camr * 1.25);
				}
			}
			
			imguiEndScrollArea();
			
		}
		
		// Test cases
		if (showTestCases)
		{
			static int testScroll = 0;
			if (imguiBeginScrollArea("Choose Test To Run", width-10-250-10-200, height-10-450, 200, 450, &testScroll))
				mouseOverMenu = true;

			vector<string>::const_iterator fileIter = files.begin();
			vector<string>::const_iterator filesEnd = files.end();
			vector<string>::const_iterator testToLoad = filesEnd;
			for (; fileIter != filesEnd; ++fileIter)
			{
				if (imguiItem(fileIter->c_str()))
				{
					testToLoad = fileIter;
				}
			}
			
			if (testToLoad != filesEnd)
			{
				string path = testCasesFolder + "/" + *testToLoad;
				test = new TestCase;
				if (test)
				{
					// Load the test.
					if (!test->load(path))
					{
						delete test;
						test = 0;
					}

					// Create sample
					Sample* newSample = 0;
					for (int i = 0; i < g_nsamples; ++i)
					{
						if (g_samples[i].name == test->getSampleName())
						{
							newSample = g_samples[i].create();
							if (newSample)
								sampleName = g_samples[i].name;
						}
					}

					delete sample;
					sample = newSample;

					if (sample)
					{
						sample->setContext(&ctx);
						showSample = false;
					}

					// Load geom.
					meshName = test->getGeomFileName();
					
					
					path = meshesFolder + "/" + meshName;
					
					delete geom;
					geom = new InputGeom;
					if (!geom || !geom->load(&ctx, path))
					{
						delete geom;
						geom = 0;
						delete sample;
						sample = 0;
						showLog = true;
						logScroll = 0;
						ctx.dumpLog("Geom load log %s:", meshName.c_str());
					}
					if (sample && geom)
					{
						sample->handleMeshChanged(geom);
					}

					// This will ensure that tile & poly bits are updated in tiled sample.
					if (sample)
						sample->handleSettings();

					ctx.resetLog();
					if (sample && !sample->handleBuild())
					{
						ctx.dumpLog("Build log %s:", meshName.c_str());
					}
					
					if (geom || sample)
					{
						const double* bmin = 0;
						const double* bmax = 0;
						if (geom)
						{
							bmin = geom->getNavMeshBoundsMin();
							bmax = geom->getNavMeshBoundsMax();
						}
						// Reset camera and fog to match the mesh bounds.
						if (bmin && bmax)
						{
							camr = sqrt(rcSqr(bmax[0] - bmin[0]) +
										 rcSqr(bmax[1] - bmin[1]) +
										 rcSqr(bmax[2] - bmin[2])) / 2;
							cameraPos[0] = (bmax[0] + bmin[0]) / 2 + camr;
							cameraPos[1] = (bmax[1] + bmin[1]) / 2 + camr;
							cameraPos[2] = (bmax[2] + bmin[2]) / 2 + camr;
							camr *= 3;
						}
						cameraEulers[0] = 45;
						cameraEulers[1] = -45;
						glFogf(GL_FOG_START, camr * 0.2);
						glFogf(GL_FOG_END, camr * 1.25);
					}
					
					// Do the tests.
					if (sample)
						test->doTests(sample->getNavMesh(), sample->getNavMeshQuery());
				}
			}				
				
			imguiEndScrollArea();
		}

		
		// Log
		if (showLog && showMenu)
		{
			if (imguiBeginScrollArea("Log", 250 + 20, 10, width - 300 - 250, 200, &logScroll))
				mouseOverMenu = true;
			for (int i = 0; i < ctx.getLogCount(); ++i)
				imguiLabel(ctx.getLogText(i));
			imguiEndScrollArea();
		}
		
		// Left column tools menu
		if (!showTestCases && showTools && showMenu) // && geom && sample)
		{
			if (imguiBeginScrollArea("Tools", 10, 10, 250, height - 20, &toolsScroll))
				mouseOverMenu = true;

			if (sample)
				sample->handleTools();
			
			imguiEndScrollArea();
		}
		
		// Marker
		if (markerPositionSet && gluProject((GLdouble)markerPosition[0], (GLdouble)markerPosition[1], (GLdouble)markerPosition[2],
								  modelviewMatrix, projectionMatrix, viewport, &x, &y, &z))
		{
			// Draw marker circle
			glLineWidth(5.0);
			glColor4ub(240,220,0,196);
			glBegin(GL_LINE_LOOP);
			const double r = 25.0;
			for (int i = 0; i < 20; ++i)
			{
				const double a = (double)i / 20.0 * RC_PI*2;
				const double fx = (double)x + cos(a)*r;
				const double fy = (double)y + sin(a)*r;
				glVertex2f(fx,fy);
			}
			glEnd();
			glLineWidth(1.0);
		}
		
		imguiEndFrame();
		imguiRenderGLDraw();		
		
		glEnable(GL_DEPTH_TEST);
		SDL_GL_SwapWindow(window);
	}
	
	imguiRenderGLDestroy();
	
	SDL_Quit();
	
	delete sample;
	delete geom;
	
	return 0;
}
