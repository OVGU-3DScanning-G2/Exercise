#ifndef MY_SHADER_H
#define MY_SHADER_H

#include <string>
#include "windows.h"
#include "gl/gl.h"
#include "glext.h"

//Shader programs are written in GLSL. GLSL is a separate programming language for graphics cards.
//The source code is provided as a text string. This string is transferred to the graphics card (driver)
//and then it becomes compiled and linked to an executable like normal C/C++-programs. But the binary code
//is not executed with CPU/RAM but with GPU/VRAM directly.

//Note: the following text strings contain the sources of the Vertex and Fragment Shader. The text string
//are enclosed in R(" text "), we as the R("") indicates a so called "string laterals" (C++11 feature), which means that 
//the compiler handles the entire text as a string even if it contains multiple lines. Without R("") each line must be
//finished with "\n" to indicate a new line.

//The programs are written to handle a single Vertex / Pixel. The GPU calls these programs in many parallel threads for
//all Vertices/Pixels.

//Here is GLSL source code for the vertex shader
static std::string VertexShaderSource = R"(
//we pass the following value from VertexShader to FragmentShader
varying vec3 N;           //varying data types can be passed from vertex shader to fragment shader. N will store the normal vector.
varying vec3 v;           //vec3, vec4 are data types integrated in GLSL, v stores later the 3d position that we define by calling glVertex or glVertexPointer, ...
varying vec4 FrontColor;  //will store the color that we define by calling glColor(r,g,b)

void main(void)                               //all the GLSL-functions that start with gl_xxxx are GLSL-functions which mainly serve as an GPU interface to our CPU-code
{
  v = vec3(gl_ModelViewMatrix * gl_Vertex);   //"gl_ModelViewMatrix" is a GLSL variable, that contains the current ModelView-Matrix (changed outside by our calls to glRotate, glTranslate, ...)
                                              //"gl_Vertex" is the 3D coordinat0e that the GPU should process (provided by our calls to glVertex3(x,y,z) or glVertexPointer(....))
  N = normalize(gl_NormalMatrix * gl_Normal); //"gl_Normal" is the 3D normal vector that the GPU should process (provided by our calls to glNormal3(x,y,z) or glNormalPointer(....))
  FrontColor = gl_Color;                      //"gl_Color" contains the color that we have set per point with glColor3ub(R,G,B) or glColorPointer.
  
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex; //looks as a duplicate, but this function must always be supplied, because the Vertex shader always uses the gl_Position variable which cannot be empty!
}
)";

//Here is GLSL source code for the fragment shader (aka "pixel shader")
//Fragment shader is for manipulating the 2D-pixels after rasterization (2d screen position and color). 
//Here we manipulate the color by applying a Phong Shading, which let our pointclouds look nice. 
//Phong shading simulates real light behaviour, where the final color results from the sum of the ambient,diffuse und specular reflection terms.
//The amount of light reflected depends on the position of camera and light source and the orientation of the surface (normal vectors). 
static std::string FragmentShaderSource = R"(
varying vec3 N;
varying vec3 v;
varying vec4 FrontColor;

void main(void)
{
  //our own definitions
  vec3 lightPosition=vec3(1.0,1.0,1.0);         //in normalized coordinates (1,1,1). I want the light coming diagonal from the back
  vec4 ambientColor=vec4(0.1, 0.1, 0.1, 1.0);   //ambient (surrounding) color is quite dark 
  //vec4 diffuseColor=vec4(0.4, 0.4, 0.4, 1.0); //you might give the diffuse light a certain fixed color, but I want to take the color that we define outside by glColor (transferred from Vertex shader using "varying" variable FrontColor)
  vec4 specularColor=vec4(0.7, 0.7, 0.7, 1.0);  //Shiny surfaces create reflecting spots, the color of these spots should be bright
  float shininess=100.0;                        //Shininess - the higher the reflective the surface (is a parameter of the specular term)
  
  vec3 N = normalize(N);                        //if not already done outside, we normalize the normal vector here

  vec3 L = normalize(lightPosition - v);        //L is the vector from the 3D position to the light position (v is 3D and was transferred from the vertex shader using a "varying" variable)
  //vec3 L = normalize( vec3(1,1,1) );

  vec3 E = normalize(-v);                       // is the vector from 3d position to the camera. We are in Eye Coordinates, so EyePos is (0,0,0) 
  vec3 R = normalize(-reflect(L, N));           // For a given incident vector I and surface normal N reflect returns the reflection direction calculated as I - 2.0 * dot(N, I) * N. N should be normalized.

  //calculate Ambient Term: 
  vec4 Iamb = ambientColor;                     // we just assign our own definition. The variable Iamb is useless it is just for sticking to the names given in the Phong formula  

  //calculate Diffuse Term:                     // the amount of diffuse reflected light depends on the surface normal N and the light position L
  vec4 Idiff = FrontColor * max(abs(dot(N, L)), 0.0);  //we take our object color as diffuse color  !! ABS was added by our own as a trick to overcome inconsistent normal orientation !!
  //vec4 Idiff=vec4(abs(N.x), abs(N.y), abs(N.z), 1)* max(abs(dot(N, L)), 0.0); 


  Idiff = clamp(Idiff, 0.0, 1.0);               //make sure that Idiff ranges between 0 and 1

  // calculate Specular Term: taken from Phong formula
  vec4 Ispec = specularColor * pow(max(dot(R, E), 0.0), 0.03* shininess);
  Ispec = clamp(Ispec, 0.0, 1.0);

  gl_FragColor = Iamb + Idiff + Ispec; //gl_FragColor is the color that the fragment (aka pixel) on the screen is assigned to. In the Phong model its the sum of ambient,diffuse and specular reflection terms
}

)";

///////////////////////////////////////////////////////////////////////////////////////////////////////
// if you are using QT 5 use existing shader classes (#include <QtOpenGL\QGLShaderProgram>)          //
// The relevant classes are QGLShader and QGLShaderProgram (see example in QT doc for these classes) //
///////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
//If you are not using QT, we have to load OpenGL-functions our own (all thats > OpenGL 1.1) //
//Therefore we implement our own shader class to simplify this process                       //
///////////////////////////////////////////////////////////////////////////////////////////////
/** @brief class that implements an easy interface to OpenGL shader functionality.

    @par usage example: 
    @code{.cpp}
      Shader myShaderProgram;                         //create only once
      myShaderProgram.initializeOpenGLextensions();   //call only once but after you have a valid OpenGL context (e.g. after the call to glfwMakeContextCurrent(window);)

      //in your render routine
      myShader.bind(); //activate
        //...draw your vertices (glVertex, glDrawArrays,....)
      myShader.unbind(); //deactivate shader

      //...draw other things that should not (or can not) be processed by this shader
    @endcode 
*/
class Shader
{
  public:
    Shader();

    void createShaderProgram(); ///< creates our shader program (if you dont have QT)
    void bind()   const { if (m_programID != 0) glUseProgram(m_programID); }  ///< activates shader program on the GPU
    void unbind() const { if (m_programID != 0) glUseProgram(0); }            ///< deactivates shader program on the GPU

  private:
    void initializeOpenGLextensions();
    bool checkShaderCompilationStatus(GLhandleARB shaderObject) const;
    GLuint m_programID = 0;

    //These cryptic types are pointers (defined in glext.h). They refer to the entry points in the openGL library DLL (graphics driver).
    PFNGLCREATEPROGRAMPROC    glCreateProgram = 0;
    PFNGLCREATESHADERPROC     glCreateShader = 0;
    PFNGLSHADERSOURCEPROC     glShaderSource = 0;
    PFNGLCOMPILESHADERPROC    glCompileShader = 0;
    PFNGLATTACHSHADERPROC     glAttachShader = 0;
    PFNGLLINKPROGRAMPROC      glLinkProgram = 0;
    PFNGLUSEPROGRAMPROC       glUseProgram = 0;
    PFNGLGETSHADERIVPROC      glGetShaderiv = 0;
    PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog = 0;
    PFNGLGETOBJECTPARAMETERIVARBPROC glGetObjectParameterivARB = 0;
};
/**/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////





#endif