//C++ RayTracer using PPM files as output

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

const int WIDTH = 800;
const int HEIGHT = 800;
const double VIEWPORT_SIZE = 1;
const double PROJECTION_PLANE_D = 1;
const double INF = numeric_limits<double>::infinity();

//Struct para colores
struct Color {
    int r, g, b;

    Color operator*(double scalar) const {
        return {
            static_cast<int>(min(255.0, r * scalar)),
            static_cast<int>(min(255.0, g * scalar)),
            static_cast<int>(min(255.0, b * scalar))
        };
    }

    Color operator+(const Color& other) const {
        return {
            min(255, r + other.r),
            min(255, g + other.g),
            min(255, b + other.b)
        };
    }
};

//Color del background de la escena (sera negro en este caso)
const Color BACKGROUND_COLOR = {0, 0, 0};

//Clase que representa un vector de 3 dimensiones (con )
class Vector3D {
public:
    double x, y, z;

    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    Vector3D() : x(0), y(0), z(0) {}

    Vector3D operator+(const Vector3D& v) const { return Vector3D(x + v.x, y + v.y, z + v.z); }//Anade 2 vectores
    Vector3D operator-(const Vector3D& v) const { return Vector3D(x - v.x, y - v.y, z - v.z); }//Resta 2 vectores
    Vector3D operator*(double scalar) const { return Vector3D(x * scalar, y * scalar, z * scalar); } //Aplica Multiplicacion escalar a los componentes

    double dot(const Vector3D& v) const { return x * v.x + y * v.y + z * v.z; } //Calcula Producto punto
    double magnitude() const { return sqrt(x * x + y * y + z * z); } //Calcula Magnitud
    Vector3D normalize() const { double mag = magnitude(); return Vector3D(x / mag, y / mag, z / mag); } //Calcula Normalisacion
};

//Structo que representa una esfera y multiples propiedades en el espacio 3D
struct Sphere {
    Vector3D center;
    double radius;
    Color color;
    double specular; // Specular exponent
    double reflective; // Reflectivity (0.0 to 1.0)

    //Constructor que inicializa las propiedades
    Sphere(Vector3D c, double r, Color col, double s, double refl) 
        : center(c), radius(r), color(col), specular(s), reflective(refl) {}
};

//Structo que representa un plano en el espacio 3D
struct Plane {
    Vector3D point;
    Vector3D normal; //vector normal del plano
    Color color;
    double specular; 
    double reflective;

    //Constructor
    Plane(Vector3D p, Vector3D n, Color col, double s, double refl)
        : point(p), normal(n.normalize()), color(col), specular(s), reflective(refl) {}
};

//Structo que representa las diferentes luces que pueden representarse en la escena
struct Light {
    string type;
    double intensity;
    Vector3D position;  
    Vector3D direction; //Para la luz direcional 
};

//Structo que ayuda a representar la posicion de la camara en la escena
struct Camera {
    Vector3D position;
    Vector3D forward; 
    Vector3D right;   
    Vector3D up;      
};

//Esferas en la escena
vector<Sphere> scene_spheres = {
    Sphere(Vector3D(2, -1.2, 3), 1, {100, 0, 255}, 500, 0.2),   //Esfera violeta
    Sphere(Vector3D(-2, -1.2, 1), 1, {0, 0, 255}, 500, 0.3),    //Esfera azul
};

//Planos en la escena
vector<Plane> scene_planes = {
    Plane(Vector3D(3, -2, 1), Vector3D(0, 1, 0), {128, 128, 128}, 500, 0.1), //Plano utilizado como superficie
    Plane(Vector3D(3, 0, 3), Vector3D(-0.5, 0, 1).normalize(), {200, 200, 200}, 300, 0.2) //Plano vertical, normalizado para que quede angulado
};

//Luces en escena
vector<Light> lights = {
    {"ambient", 0.2, Vector3D(), Vector3D()},
    {"point", 0.6, Vector3D(2, 1, 0), Vector3D()},
    {"directional", 0.2, Vector3D(), Vector3D(1, 4, 4).normalize()}
};

//Posicion de la camara en la escena
Camera camera = {
    Vector3D(0, 0, -5),
    Vector3D(0, 0, 1), 
    Vector3D(1, 0, 0), 
    Vector3D(0, 1, 0) 
};

//Funcion CanvasToViewport que convierte cordenadas de pixeles del canvas a coordenas del viewport en el espacio 3D
Vector3D CanvasToViewport(int x, int y) {
    return Vector3D(x * VIEWPORT_SIZE / WIDTH, y * VIEWPORT_SIZE / HEIGHT, PROJECTION_PLANE_D);
}

//Funcion donde transforma la direcion en el espacio del viewport a un espacio en la escena como tal
Vector3D ApplyCameraRotation(const Vector3D& v) {
   
    return Vector3D( //Aplica rotacion de camara al vector
        v.x * camera.right.x + v.y * camera.up.x + v.z * camera.forward.x,
        v.x * camera.right.y + v.y * camera.up.y + v.z * camera.forward.y,
        v.x * camera.right.z + v.y * camera.up.z + v.z * camera.forward.z
    );
}

//Funcion que calcula donde el rayo interseca con una esfera en el espacio 3D
pair<double, double> IntersectRaySphere(Vector3D O, Vector3D D, Sphere sphere) {
    Vector3D CO = O - sphere.center;

    double a = D.dot(D);
    double b = 2 * CO.dot(D);
    double c = CO.dot(CO) - sphere.radius * sphere.radius;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return {INF, INF};

    double t1 = (-b + sqrt(discriminant)) / (2 * a);
    double t2 = (-b - sqrt(discriminant)) / (2 * a);

    return {t1, t2};
}

//Funcion quye calcula donde el rayo interseca con un plano en el espacio 3D
double IntersectRayPlane(Vector3D O, Vector3D D, Plane plane) {
    double denom = plane.normal.dot(D);
    if (abs(denom) < 1e-6) {
        return INF; // Ray is parallel to the plane
    }
    Vector3D PO = plane.point - O;
    double t = PO.dot(plane.normal) / denom;
    return (t >= 0) ? t : INF; // Return t only if it's positive
}

//Se busca el objeto mas cercano que un rayo interseca en el espacio 3D
pair<void*, double> ClosestIntersection(Vector3D O, Vector3D D, double t_min, double t_max) {
    double closest_t = INF;
    void* closest_object = nullptr;
    bool is_sphere = true;

    //Verificar intersecciones con esferas
    for (auto& sphere : scene_spheres) {
        auto [t1, t2] = IntersectRaySphere(O, D, sphere);

        if (t1 >= t_min && t1 <= t_max && t1 < closest_t) {
            closest_t = t1;
            closest_object = &sphere;
            is_sphere = true;
        }
        if (t2 >= t_min && t2 <= t_max && t2 < closest_t) {
            closest_t = t2;
            closest_object = &sphere;
            is_sphere = true;
        }
    }

    //Verificar intersecciones con planos
    for (auto& plane : scene_planes) {
        double t = IntersectRayPlane(O, D, plane);
        if (t >= t_min && t <= t_max && t < closest_t) {
            closest_t = t;
            closest_object = &plane;
            is_sphere = false;
        }
    }

    return {closest_object, closest_t};
}

//Funcion que calcula la direcfcion de la reflecion de un rayo que brinca con una superficie
Vector3D ReflectRay(Vector3D R, Vector3D N) {
    return N * 2.0 * N.dot(R) - R;
}

//Funcion que calcula la intensidad de la luz en un punto en la escena
double ComputeLighting(Vector3D P, Vector3D N, Vector3D V, double s) {
    double i = 0.0;
    for (const auto& light : lights) {
        if (light.type == "ambient") {
            i += light.intensity;
        } else {
            Vector3D L;
            double t_max;
            if (light.type == "point") {
                L = light.position - P;
                t_max = 1;
            } else {
                L = light.direction;
                t_max = INF;
            }

            
            auto [shadow_object, shadow_t] = ClosestIntersection(P, L, 0.001, t_max);
            if (shadow_object != nullptr) {
                continue; 
            }

            
            double n_dot_l = N.dot(L);
            if (n_dot_l > 0) {
                i += light.intensity * n_dot_l / (N.magnitude() * L.magnitude());
            }

            
            if (s != -1) {
                Vector3D R = ReflectRay(L, N); 
                double r_dot_v = R.dot(V);
                if (r_dot_v > 0) {
                    i += light.intensity * pow(r_dot_v / (R.magnitude() * V.magnitude()), s);
                }
            }
        }
    }
    return i;
}

//Funcion que representa un rayo que traza en la escena y que determina que color debe ser representado para cada pixel
Color TraceRay(Vector3D O, Vector3D D, double t_min, double t_max, int recursion_depth) {
    auto [closest_object, closest_t] = ClosestIntersection(O, D, t_min, t_max);

    if (closest_object == nullptr) {
        return BACKGROUND_COLOR;
    }

    //Computar punto de interseccion
    Vector3D P = O + D * closest_t;

    //Computar el color y la normal dependiendo el tipo de objeto
    Vector3D N;
    Color color;
    double specular, reflective;

    //Determinar si el objeto mas cerca es una esfera o plano
    bool is_sphere = false;
    for (auto& sphere : scene_spheres) {
        if (&sphere == closest_object) {
            is_sphere = true;
            N = (P - sphere.center).normalize();
            color = sphere.color;
            specular = sphere.specular;
            reflective = sphere.reflective;
            break;
        }
    }

    if (!is_sphere) {
        for (auto& plane : scene_planes) {
            if (&plane == closest_object) {
                N = plane.normal;
                color = plane.color;
                specular = plane.specular;
                reflective = plane.reflective;
                break;
            }
        }
    }

    //Computar la luz
    Vector3D V = D * -1.0; 
    double intensity = ComputeLighting(P, N, V, specular);
    Color local_color = color * intensity;

    //Limite de recursion
    if (recursion_depth <= 0 || reflective <= 0) {
        return local_color;
    }

    //Computar color reflejado
    Vector3D R = ReflectRay(D * -1.0, N);
    Color reflected_color = TraceRay(P, R, 0.001, INF, recursion_depth - 1);

    //Mezclar color local y color reflejado (basado en la reflectividad)
    return local_color * (1 - reflective) + reflected_color * reflective;
}

//Funcion Main
int main() {
    ofstream image("output.ppm");
    image << "P3\n" << WIDTH << " " << HEIGHT << "\n255\n";

    for (int y = -HEIGHT / 2; y < HEIGHT / 2; ++y) {
        for (int x = -WIDTH / 2; x < WIDTH / 2; ++x) {
            // Compute ray direction in viewport space
            Vector3D D_viewport = CanvasToViewport(x, -y);

            // Apply camera rotation to the ray direction
            Vector3D D = ApplyCameraRotation(D_viewport).normalize();

            // Trace the ray from the camera's position
            Color color = TraceRay(camera.position, D, 1, INF, 3); // Recursion depth of 3
            image << color.r << " " << color.g << " " << color.b << " ";
        }
    }

    image.close();

    return 0;
}