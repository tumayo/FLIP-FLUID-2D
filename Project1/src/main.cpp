#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "gluvi.h"
#include "fluidsim.h"
#include "openglutils.h"
#include "array2_utils.h"

#include "fluidsim.cpp"
#include "gluvi.cpp"
#include "openglutils.cpp"

//Comment this to disable frame dumping
#define MAKE_MOVIE


using namespace std;


//Try changing the grid resolution
int grid_resolution = 30;
float timestep = 0.001f;
int pause_btw_frames_in_ms = 10;
char img_file_path[]{ "C:/output/" };

//Display properties
bool draw_grid = true;
bool draw_particles = true;
bool draw_velocities = true;
bool draw_boundaries = true;

float grid_width = 1;

FluidSim sim;

//Gluvi stuff
//-------------
Gluvi::PanZoom2D cam(-0.1f, -0.35f, 1.2f);
double oldmousetime;
Vec2f oldmouse;
void display();
void mouse(int button, int state, int x, int y);
void drag(int x, int y);
void timer(int junk);


//Boundary definition - several circles in a circular domain.

Vec2f c0(0.5f, 0.5f), c1(0.7f, 0.5f), c2(0.3f, 0.35f), c3(0.5f, 0.7f);
float rad0 = 0.4f, rad1 = 0.1f, rad2 = 0.1f, rad3 = 0.1f;

float circle_phi(const Vec2f& position, const Vec2f& centre, float radius) {
    return (dist(position, centre) - radius);
}

float boundary_phi(const Vec2f& position) {
    float phi0 = -circle_phi(position, c0, rad0);
    //float phi1 = circle_phi(position, c1, rad1);
    //float phi2 = circle_phi(position, c2, rad2);
    //float phi3 = circle_phi(position, c3, rad3);
    //return min(min(phi0,phi1),min(phi2,phi3));
    return phi0;
}

//Draw velocity for each particle -- Tumay
void draw_velocity() {
    for (int i = 0; i < sim.particles.size(); ++i) {
        Vec2f pos = sim.particles[i];
        Vec2f vel = sim.get_velocity(sim.particles[i]);

        float vel_norm = sqrt(pow(vel[0], 2) + pow(vel[1], 2));
        Vec2f end = Vec2f(pos[0] + vel[0] / vel_norm / 20, pos[1] + vel[1] / vel_norm / 20);

        glBegin(GL_POINTS);
        glColor3f(1, 1, 1);
        glVertex2f(end[0], end[1]);
        glEnd();

        glColor3f(1, 0, 0);
        glBegin(GL_LINES);
        glVertex2f(pos[0], pos[1]);
        glVertex2f(end[0], end[1]);
        glEnd();
    }
}



//Main testing code
//-------------
int main(int argc, char** argv)
{

    //Setup viewer stuff
    Gluvi::init("Basic Fluid Solver with Vorticity", &argc, argv);
    Gluvi::camera = &cam;
    Gluvi::userDisplayFunc = display;
    Gluvi::userMouseFunc = mouse;
    Gluvi::userDragFunc = drag;
    glClearColor(1, 1, 1, 1);

    glutTimerFunc(1000, timer, 0);

    //Set up the simulation
    sim.initialize(grid_width, grid_resolution, grid_resolution);
    sim.set_boundary(boundary_phi);

    /*for(int i = 0; i < sqr(grid_resolution); ++i) {
       float x = randhashf(i*2, 0,1);
       float y = randhashf(i*2+1, 0,1);
       Vec2f pt(x,y);
       if(boundary_phi(pt) > 0)
          sim.add_particle(pt);
    }*/
    float dx = grid_width / (float)grid_resolution;
    for (int i = 0; i < grid_resolution + 1; i++) {
        for (int j = 0; j < grid_resolution + 1; j++) {
            float x = i * dx;
            float y = j * dx;
            if (0.4 <= x && x <= 0.7 && 0.4 <= y && y <= 0.7) {
                sim.add_particle(Vec2f(x, y));
            }
        }
    }
    //sim.advance(timestep);
    Gluvi::run();

    return 0;
}


void display(void)
{

    if (draw_grid) {
        glColor3f(0, 0, 0);
        glLineWidth(1);
        //draw_grid2d(Vec2f(0,0), sim.dx, sim.ni, sim.nj);  
    }

    if (draw_boundaries) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        draw_circle2d(c0, rad0, 50);
        /*draw_circle2d(c1, rad1, 50);
        draw_circle2d(c2, rad2, 50);
        draw_circle2d(c3, rad3, 50);

        //There's a bug, so draw one more(?)
        draw_circle2d(c3, 0, 10);*/
    }

    if (draw_particles) {
        glColor3f(0, 0, 0);
        glPointSize(3);
        draw_points2d(sim.particles);
    }

    /*if(draw_velocities) {
       for(int j = 0;j < sim.nj; ++j) for(int i = 0; i < sim.ni; ++i) {
          Vec2f pos((i+0.5f)*sim.dx,(j+0.5f)*sim.dx);
          draw_arrow2d(pos, pos + 0.01f*sim.get_velocity(pos), 0.1f*sim.dx);
          //cout << "Vel x: " << sim.get_velocity(pos)[0] << "Vel y:" << sim.get_velocity(pos)[1] << endl;
       }
    }*/

    draw_velocity();
    /*glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int i = 0; i <= sim.ni; i++) {
        for (int j = 0; j < sim.nj; j++) {

            if (sim.vor(i, j) > 10)
                glColor3f(1, 0, 0);
            else if (sim.vor(i, j) > 5 && sim.vor(i, j) < 10)
                glColor3f(0.66, 0, 0);
            else if (sim.vor(i, j) < 5 && sim.vor(i,j) > 0)
                glColor3f(0.33, 0, 0);
            else if (sim.vor(i, j) < -10)
                glColor3f(0, 0, 1);
            else if (sim.vor(i, j) > -10 && sim.vor(i, j) < -5)
                glColor3f(0, 0, 0.66);
            else if (sim.vor(i, j) < 0 && sim.vor(i, j) > -5)
                glColor3f(0, 0, 0.33);
            else
                glColor3f(0, 1, 0);
            glVertex2f(i * sim.dx, j * sim.dx);
        }
    }
    glEnd();*/

}

void mouse(int button, int state, int x, int y)
{
    Vec2f newmouse;
    cam.transform_mouse(x, y, newmouse.v);
    //double newmousetime=get_time_in_seconds();

    oldmouse = newmouse;
    //oldmousetime=newmousetime;
}

void drag(int x, int y)
{
    Vec2f newmouse;
    cam.transform_mouse(x, y, newmouse.v);
    //double newmousetime=get_time_in_seconds();

    oldmouse = newmouse;
    //oldmousetime=newmousetime;
}

void timer(int junk)
{
    sim.advance(timestep);


#ifdef MAKE_MOVIE
    static int frame = 0;
    frame++;

    char* sgifileformat;
    sgifileformat = new char[strlen(img_file_path) + 50];

    sprintf(sgifileformat, "%s/screenshot%%04d.sgi", img_file_path);
    Gluvi::sgi_screenshot(sgifileformat, frame);

    delete[] sgifileformat;
#endif

    glutPostRedisplay();
    glutTimerFunc(pause_btw_frames_in_ms, timer, 0);

}





