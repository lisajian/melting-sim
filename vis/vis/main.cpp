#define OUTPUT_ANIMATION 1

#include <stdlib.h>
#include <stdio.h>

#include "Particles.h"
#include "FluidSim.h"
#include "FluidSim.cpp"

#if OUTPUT_ANIMATION
#include <opencv2/opencv.hpp>
#endif


inline float clip(const float& n, const float& lower, const float& upper)
{
    return glm::max(lower, glm::min(n, upper));
}

float theta = M_PI/8;
float phi = -M_PI/8+M_PI_2;
float dist = 2.5;
int width = 800;
int height = 800;
int frame = 0;
const int render_step = 1;
int mx, my;
float _zoom = 1.0;

// Particles particles;
FluidSim fluidSim;

void display(void);

void reshape(int width, int height);

void idle(void) {
    glutPostRedisplay();
    return;
    if (frame/render_step < 10) {
        std::cout << "step: " << frame << std::endl;
        fluidSim.step();
        glutPostRedisplay();
    }
    // fluidSim.step();
    // glutPostRedisplay();
    if(frame/render_step >= 10) {
        // std::cout << "stopped" << std::endl;
        glutPostRedisplay();
        return;
    }
    if(frame%render_step == 0)
    {
        // std::cout << frame << std::endl;
        #if OUTPUT_ANIMATION
        cv::Mat3b image(height, width);
        glReadPixels(0, 0, width, height, GL_BGR, GL_UNSIGNED_BYTE, image.data);
        cv::flip(image, image, 0);
        char fn[512];
        sprintf(fn, "result/%04d.png", frame/render_step);
        cv::imwrite(fn, image);
        #endif
    }
    frame++;
}

void mouse(int button, int state, int x, int y);

void motion(int x, int y);

void keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case 'r' :
            frame = 0;
            fluidSim.reset();
            break;
        case '=':
            _zoom += 0.3;
            break;
        case '-':
            _zoom -= 0.3;
            if (_zoom <= 0.1) {
                _zoom = 0.1;
            }
            break;
        case 'n':
            fluidSim.step();
            glutPostRedisplay();
            break;
        default:
            break;
    }
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(width, height);

    (void)glutCreateWindow("GLUT Program");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutMainLoop();

    return EXIT_SUCCESS;
}

void reshape(int w, int h) {
    width = w;
    height = h;
    glViewport(0, 0, w, h);
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // your drawing code goes here
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90, 1, 0.01, 100);
    gluLookAt(dist*sin(phi)*cos(theta), dist*cos(phi), dist*sin(phi)*sin(theta),
            0, 0, 0,
            0, 1, 0);
    glScalef(_zoom,_zoom,_zoom);

    // bottom face
    glBegin(GL_LINE_LOOP);
    glVertex3f(-2, -2, -2);
    glVertex3f(-2, 2, -2);
    glVertex3f(2, 2, -2);
    glVertex3f(2, -2, -2);
    glEnd();

    // right face
    glBegin(GL_LINE_LOOP);
    glVertex3f(2, 2, -2);
    glVertex3f(2, -2, -2);
    glVertex3f(2, -2, 2);
    glVertex3f(2, 2, 2);
    glEnd();

    // front face
    glBegin(GL_LINE_LOOP);
    glVertex3f(-2, 2, -2);
    glVertex3f(2, 2, -2);
    glVertex3f(2, 2, 2);
    glVertex3f(-2, 2, 2);
    glEnd();

    // left face
    glBegin(GL_LINE_LOOP);
    glVertex3f(-2, 2, -2);
    glVertex3f(-2, -2, -2);
    glVertex3f(-2, -2, 2);
    glVertex3f(-2, 2, 2);
    glEnd();

    // back face
    glBegin(GL_LINE_LOOP);
    glVertex3f(-2, -2, -2);
    glVertex3f(2, -2, -2);
    glVertex3f(2, -2, 2);
    glVertex3f(-2, -2, 2);
    glEnd();

    // top face
    glBegin(GL_LINE_LOOP);
    glVertex3f(2, 2, 2);
    glVertex3f(2, -2, 2);
    glVertex3f(-2, -2, 2);
    glVertex3f(-2, 2, 2);
    glEnd();

    fluidSim.render();

    glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {
    if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        mx = x;
        my = y;
    }
}

void motion(int x, int y) {
    int dx = x - mx;
    int dy = y - my;
    mx = x;
    my = y;
    if(abs(dx) > abs(dy))
        theta += dx*0.005;
    else
        phi -= dy*0.005;
    if(theta > 2*M_PI)
        theta -= 2*M_PI;
    if(theta < -2*M_PI)
        theta += 2*M_PI;
    phi = clip(phi, M_PI/12, M_PI*11/12);
    glutPostRedisplay();
}
