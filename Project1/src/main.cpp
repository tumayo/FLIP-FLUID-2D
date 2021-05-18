#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>
#include <vector>

#include "gluvi.h"
#include "fluidsim.h"
#include "openglutils.h"
#include "array2_utils.h"
#include "colorf.h"

#include "fluidsim.cpp"
#include "gluvi.cpp"
#include "openglutils.cpp"

#include "makelevelset2.h"

//Comment this to disable frame dumping
#define MAKE_MOVIE

using namespace std;

//Try changing the grid resolution
int grid_resolution = 40;
float timestep = 0.05f; 
int pause_btw_frames_in_ms = 10;
char img_file_path[]{ "C:/output/" };

//Display properties
bool draw_grid = true;
bool draw_particles = true;
bool draw_velocities = true;
bool draw_boundaries = true;

float grid_width = 1;
int start = 0;

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
void key_func(unsigned char key, int x, int y);

//Tumay -- Bunny Boundary Construct
const int num_points = 60; //# of vertices in the bunny boundary
float temp_BC[num_points * 2] = {};
vector<Vec2f> BC_vertices(num_points, 0); //list of vertices in the bunny boundary
vector<Vec2ui> BC_edges(num_points, 0); //list of edges in the bunny boundary

void BC_construct(void);
void draw_BC();

//Tumay -- Box Boundary Construct
void Box_BC_construct(float dx);
vector<Vec2f> Box_BC_vertices(4, 0); 
vector<Vec2ui> Box_BC_edges(4, 0);

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

//Tumay -- Draw velocity for each particle
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

void draw_nodal_solid_phi() {
    for (int i = 0; i < sim.ni + 1; i++) {
        for (int j = 0; j < sim.nj + 1; j++) {
            glBegin(GL_POINTS);
            if (sim.nodal_solid_phi(i, j) < 0) {
                glColor3f(1, 1, 0);
                glVertex2f(i * sim.dx, j * sim.dx);
            }
            glEnd();
        }
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
    glutKeyboardFunc(key_func);

    //Set up the simulation
    sim.initialize(grid_width, grid_resolution, grid_resolution);
    //Existing Circle Boundary
    //sim.set_boundary(boundary_phi);

   /*Box_BC_construct(sim.dx);
    Vec2f origin = Vec2f(0.0, 0.0);
    make_level_set2(Box_BC_edges, Box_BC_vertices, origin, sim.dx, sim.ni + 1, sim.nj + 1, sim.nodal_solid_phi);
    for (int i = 0; i < sim.nodal_solid_phi.a.size(); ++i)
        sim.nodal_solid_phi.a[i] = -sim.nodal_solid_phi.a[i];*/
    
    //revert to no boundaries
    sim.nodal_solid_phi.assign(+1);

    //Tumay -- Bunny Boundary
    /*BC_construct();
    Vec2f origin = Vec2f(0.0, 0.0);
    make_level_set2(BC_edges, BC_vertices, origin, sim.dx, sim.ni+1 , sim.nj+1, sim.nodal_solid_phi);
    for (int i = 0; i < sim.nodal_solid_phi.a.size(); ++i)
       sim.nodal_solid_phi.a[i] = -sim.nodal_solid_phi.a[i];
    */

    /*for(int i = 0; i < sqr(grid_resolution); ++i) {
       float x = randhashf(i*2, 0,1);
       float y = randhashf(i*2+1, 0,1);
       Vec2f pt(x,y);
       if(boundary_phi(pt) > 0)
          sim.add_particle(pt);
    }*/

    float dx = grid_width / (float)grid_resolution;
    for (int i = 0; i < grid_resolution*2 + 1; i++) {
        for (int j = 0; j < grid_resolution*2 + 1; j++) {
            float x = i * dx/2;
            float y = j * dx/2;
            /*if(interpolate_value(Vec2f(x,y) / dx, sim.nodal_solid_phi) > 0)
                sim.add_particle(Vec2f(x, y));*/
            if (0.2f <= x && x <= 0.6f && 0.2f <= y && y <= 0.6f)
                sim.add_particle(Vec2f(x, y));
            
        }
    }
    Gluvi::run();

    return 0;
}


void display(void)
{
    if (draw_boundaries) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        //Circle Boundary
        //draw_circle2d(c0, rad0, 50);
        /*draw_circle2d(c1, rad1, 50);
        draw_circle2d(c2, rad2, 50);
        draw_circle2d(c3, rad3, 50);

        //There's a bug, so draw one more(?)
        draw_circle2d(c3, 0, 10);*/
    }
    
    /*if (draw_particles) {
        glColor3f(0, 0, 0);
        glPointSize(3);
        draw_points2d(sim.particles);
    }*/

    //Particles with interpolated vorticities
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    for (int i = 0; i < sim.particles.size(); ++i) {
        Vec2f particle_pos = sim.particles[i].v;
        float vor_value = interpolate_value(particle_pos / sim.dx, sim.vor);

        // side of the initial square length
        float side = 0.4f;
        // scale of the vorticity field
        float strength = 10.0f;
        // Max value of the vorticity on the square, for color normalization
        float norm_factor = strength * sqrt(2 * (side) * (side));

        // normalize the vorticity bewteen 0 and 1 ( < 0.5 now correspond to negative value)
        float vort_norm = 0.5f * (vor_value / norm_factor + 1.0f);
        // 1.0f - vort_norm because I mistakenly switch the colors in my simulation
        Colorf color = clamp(ramp(PortalW, 1.0f - vort_norm));
        // give to Opengl, alpha can be wathever you want
        glColor4f(color.r, color.g, color.b, 0.5);
        glVertex2fv(sim.particles[i].v);
    }
    glEnd();

    if (draw_grid) {
        glColor3f(0, 0, 0);
        glLineWidth(1);
        draw_grid2d(Vec2f(0, 0), sim.dx, sim.ni, sim.nj);
    }
  
    /*if(draw_velocities) {
       for(int j = 0;j < sim.nj; ++j) for(int i = 0; i < sim.ni; ++i) {
          Vec2f pos((i+0.5f)*sim.dx,(j+0.5f)*sim.dx);
          draw_arrow2d(pos, pos + 0.01f*sim.get_velocity(pos), 0.1f*sim.dx);
          //cout << "Vel x: " << sim.get_velocity(pos)[0] << "Vel y:" << sim.get_velocity(pos)[1] << endl;
       }
    }*/
   
    //draw_velocity();

    /*glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < sim.ni; i++) {
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

    //Bunny Boundary
    //draw_BC();
    //draw_nodal_solid_phi();
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

//Tumay 
void key_func(unsigned char key, int x, int y) {

    switch (key) { 
    case 's': start = !start; 
        glutTimerFunc(pause_btw_frames_in_ms, timer, 0);
        break;
    }
}

void timer(int junk)
{
    if (start) {
       std::cout << "Running\n";
        sim.advance(timestep);
        //sim.flip_adv_advance(timestep);

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

}

//Tumay -- Construct Bunny Boundary
static void BC_construct(void) {
    std::ifstream ifs("./figure_bunny.txt");
    std::string str;
    float a;
    int count = 0;
    while (getline(ifs, str)) {
        a = stof(str);
        temp_BC[count] = a;
        count += 1;
    }
    for (int i = 0; i < num_points; i++) {
        BC_vertices[i] = Vec2f(temp_BC[2 * i], temp_BC[2 * i + 1]);
        if (i != num_points - 1) {
            BC_edges[i] = Vec2ui(i, i + 1);
        }
        else {
            BC_edges[i] = Vec2ui(i, 0);
        }
    }
}

static void draw_BC(void)
{
    glLineWidth(5.0);
    glBegin(GL_LINES);
    for (int i = 0; i < num_points; i++) {
        glColor3f(1.0, 1.0, 0.0);
        glVertex2f(BC_vertices[i][0], BC_vertices[i][1]);
        glVertex2f(BC_vertices[(i + 1) % num_points][0], BC_vertices[(i + 1) % num_points][1]);
    }
    glEnd();
}

static void Box_BC_construct(float dx) {
    //Starting the boundary cells 2 grid cells inside the edge cells.
    dx = 2 * dx;
    Box_BC_vertices[0] = Vec2f(0.0f + dx, 0.0f + dx);
    Box_BC_vertices[1] = Vec2f(0.0f + dx, 1.0f - dx);
    Box_BC_vertices[2] = Vec2f(1.0f - dx, 1.0f - dx);
    Box_BC_vertices[3] = Vec2f(1.0f - dx, 0.0f + dx);

    Box_BC_edges[0] = Vec2ui(0, 1);
    Box_BC_edges[1] = Vec2ui(1, 2);
    Box_BC_edges[2] = Vec2ui(2, 3);
    Box_BC_edges[3] = Vec2ui(3, 0);

}





