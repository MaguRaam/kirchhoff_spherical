#include "../include/cell_properties.h"

// Class member function definitions 

// Default constructor 

cell_properties::cell_properties() {

    area = 0.0;
    cell_quadrature_points = new Point<2>[0];
    jxws = new double[0];              
    face_quadrature_points1 = new Point<2>[0];
    face_quadrature_points2 = new Point<2>[0];
    face_jxws_1 = new double[0];              
    face_jxws_2 = new double[0];              
    face_normal_x1 = new double[0]; 
    face_normal_x2 = new double[0];
    face_normal_y1 = new double[0];
    face_normal_y2 = new double[0];
    face_center_quadrature_points = new Point<2>[0];
    face_center_normal_x = new double[0];
    face_center_normal_y = new double[0];
    surface_area = new double[0]; 
	h_size = 0.0;
}


// Main constructor taking arguments 

cell_properties::cell_properties(double _area,
	Point<2> _cell_quadrature_points[], double* _jxws, 
	Point<2> _face_quadrature_points1[], Point<2> _face_quadrature_points2[],double* _face_jxws_1,
	double* _face_jxws_2, double* _face_normal_x1, double* _face_normal_x2, double* _face_normal_y1, 
	double* _face_normal_y2, Point<2> _face_center_quadrature_points[],double* _face_center_normal_x, 
	double* _face_center_normal_y, double* _surface_area) {

    area = _area; 
    
    // Initialize all the arrays 
    cell_quadrature_points = new Point<2>[4];
    jxws = new double[4];     
    face_quadrature_points1 = new Point<2>[4];
    face_quadrature_points2 = new Point<2>[4];
    face_jxws_1 = new double[4];     
    face_jxws_2 = new double[4];     
    face_normal_x1 = new double[4]; 
    face_normal_x2 = new double[4];
    face_normal_y1 = new double[4];
    face_normal_y2 = new double[4];
    face_center_quadrature_points = new Point<2>[4];
    face_center_normal_x = new double[4];
    face_center_normal_y = new double[4];
    surface_area = new double[4]; 
	h_size = std::sqrt(area);
    
    // Fill all the arrays 
    
    for (unsigned int f = 0; f < 4; f++) {

        cell_quadrature_points[f] = _cell_quadrature_points[f];
        jxws[f] = _jxws[f];    
        face_quadrature_points1[f] = _face_quadrature_points1[f];
        face_quadrature_points2[f] = _face_quadrature_points2[f];
        face_jxws_1[f] = _face_jxws_1[f]; 
        face_jxws_2[f] = _face_jxws_2[f];    
        face_normal_x1[f]     = _face_normal_x1[f]; 
        face_normal_x2[f]     = _face_normal_x2[f];
        face_normal_y1[f]     = _face_normal_y1[f];
        face_normal_y2[f]     = _face_normal_y2[f];
        face_center_quadrature_points[f] = _face_center_quadrature_points[f];
        face_center_normal_x[f]     = _face_center_normal_x[f];
        face_center_normal_y[f]     = _face_center_normal_y[f];
        surface_area[f]       = _surface_area[f]; 
    }

}

// Copy constructor 

cell_properties::cell_properties(const cell_properties& cp) {

    area = cp.area; 
    
        // Initialize all the arrays 
    
    cell_quadrature_points = new Point<2>[4];
    jxws = new double[4];     
    face_quadrature_points1 = new Point<2>[4];
    face_quadrature_points2 = new Point<2>[4];
    face_jxws_1 = new double[4];     
    face_jxws_2 = new double[4];
    face_normal_x1 = new double[4]; 
    face_normal_x2 = new double[4];
    face_normal_y1 = new double[4];
    face_normal_y2 = new double[4];
    face_center_quadrature_points = new Point<2>[4];
    face_center_normal_x = new double[4];
    face_center_normal_y = new double[4];
    surface_area = new double[4]; 
	
	h_size = std::sqrt(area);
    
    // Fill all the arrays 
    
    for (unsigned int f = 0; f < 4; f++) {
    
    	cell_quadrature_points[f] = cp.cell_quadrature_points[f];
	    jxws[f] = cp.jxws[f];      
        face_quadrature_points1[f] = cp.face_quadrature_points1[f];        
        face_quadrature_points2[f] = cp.face_quadrature_points2[f];
        face_jxws_1[f] = cp.face_jxws_1[f]; 
        face_jxws_2[f] = cp.face_jxws_2[f]; 
        face_normal_x1[f]     = cp.face_normal_x1[f]; 
        face_normal_x2[f]     = cp.face_normal_x2[f];
        face_normal_y1[f]     = cp.face_normal_y1[f];
        face_normal_y2[f]     = cp.face_normal_y2[f];
        face_center_quadrature_points[f] = cp.face_center_quadrature_points[f];
        face_center_normal_x[f]     = cp.face_center_normal_x[f];
        face_center_normal_y[f]     = cp.face_center_normal_y[f];
        surface_area[f]       = cp.surface_area[f]; 
    }
    
    
}

// Destructor 

cell_properties::~cell_properties() {

    delete[] cell_quadrature_points;
    delete[] jxws;     
    delete[] face_quadrature_points1;
    delete[] face_quadrature_points2;
    delete[] face_jxws_1;     
    delete[] face_jxws_2;     
    delete[] face_normal_x1; 
    delete[] face_normal_x2;
    delete[] face_normal_y1;
    delete[] face_normal_y2;
    delete[] face_center_quadrature_points;
    delete[] face_center_normal_x;
    delete[] face_center_normal_y;
    delete[] surface_area; 
}

// Reinitialize the object 

void cell_properties::reinit(double _area,
	 Point<2>* _cell_quadrature_points, double* _jxws, 
	 Point<2>* _face_quadrature_points1, Point<2>* _face_quadrature_points2, double* _face_jxws_1, 
	double* _face_jxws_2, double* _face_normal_x1, double* _face_normal_x2, double* _face_normal_y1, 
	double* _face_normal_y2, Point<2> _face_center_quadrature_points[],double* _face_center_normal_x, 
	double* _face_center_normal_y,double* _surface_area) {
    
    // Delete the old data 

    delete[] cell_quadrature_points;
    delete[] jxws;     
    delete[] face_quadrature_points1;
    delete[] face_quadrature_points2;
    delete[] face_jxws_1;     
    delete[] face_jxws_2;
    delete[] face_normal_x1; 
    delete[] face_normal_x2;
    delete[] face_normal_y1;
    delete[] face_normal_y2;
    delete[] face_center_quadrature_points;
    delete[] face_center_normal_x;
    delete[] face_center_normal_y;
    delete[] surface_area; 
    
    // Create the data again 
    
    area = _area; 

    // Initialize all the arrays 
    
    cell_quadrature_points = new Point<2>[4];
    jxws = new double[4];     
    face_quadrature_points1 = new Point<2>[4];
    face_quadrature_points2 = new Point<2>[4];
    face_jxws_1 = new double[4];     
    face_jxws_2 = new double[4];
    face_normal_x1 = new double[4]; 
    face_normal_x2 = new double[4];
    face_normal_y1 = new double[4];
    face_normal_y2 = new double[4];
    face_center_quadrature_points = new Point<2>[4];
    face_center_normal_x = new double[4];
    face_center_normal_y = new double[4];
    surface_area = new double[4]; 
	
	h_size = std::sqrt(area);

    // Fill all the arrays 
    
    for (unsigned int f = 0; f < 4; f++) {
    
        cell_quadrature_points[f] = _cell_quadrature_points[f];
        jxws[f] = _jxws[f];    
        face_quadrature_points1[f] = _face_quadrature_points1[f];
        face_quadrature_points2[f] = _face_quadrature_points2[f];
        face_jxws_1[f] = _face_jxws_1[f]; 
        face_jxws_2[f] = _face_jxws_2[f];
        face_normal_x1[f]     = _face_normal_x1[f]; 
        face_normal_x2[f]     = _face_normal_x2[f];
        face_normal_y1[f]     = _face_normal_y1[f];
        face_normal_y2[f]     = _face_normal_y2[f];
        face_center_quadrature_points[f] = _face_center_quadrature_points[f];
        face_center_normal_x[f]     = _face_center_normal_x[f];
        face_center_normal_y[f]     = _face_center_normal_y[f];
        surface_area[f]       = _surface_area[f]; 
    }

}


// Get the propeties of the cell 

// 1. area of the cell 

double cell_properties::measure() const {
    return area;
}

// 2. cell quadrature point

Point<2> cell_properties::cell_quadrature_point(unsigned int quad_point) const {

    assert(quad_point < 4);
    return cell_quadrature_points[quad_point]; 
}

// 3. cell jxw at quadrature point

double cell_properties::jxw(unsigned int quad_point) const {

    assert(quad_point < 4);
    return jxws[quad_point]; 
}

// 4. First quadrature point on face f 

Point<2> cell_properties::face_quadrature_point1(unsigned int face_index) const {

    assert(face_index < 4);
    return face_quadrature_points1[face_index]; 
}

// 5. Second quadrature point on face f 

Point<2> cell_properties::face_quadrature_point2(unsigned int face_index) const {

    assert(face_index < 4);
    return face_quadrature_points2[face_index]; 
}

// 6. jxw at first quadrature point on face f 

double cell_properties::face_jxw_1(unsigned int face_index) const {

    assert(face_index < 4);
    return face_jxws_1[face_index]; 
}

// 7. jxw at second quadrature point on face f 

double cell_properties::face_jxw_2(unsigned int face_index) const {

    assert(face_index < 4);
    return face_jxws_2[face_index]; 
}


// 8. First normal (x) on face 1 

double cell_properties::nx1(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_x1[face_index]; 
}

// 9. Second normal (x) on face 2 

double cell_properties::nx2(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_x2[face_index]; 
}

// 10. First normal (y) on face 1 

double cell_properties::ny1(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_y1[face_index]; 
}

// 11. Second normal (y) on face 2 

double cell_properties::ny2(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_y2[face_index]; 
}

// 12. quadrature point at center of face f 

Point<2> cell_properties::face_center_quadrature_point(unsigned int face_index) const {

    assert(face_index < 4);
    return face_center_quadrature_points[face_index]; 
}

// 13. x componenet of normal at center of face f 

double cell_properties::center_nx(unsigned int face_index) const {

    assert(face_index < 4);
    return face_center_normal_x[face_index]; 
}

// 14. y componenet of normal at center of face f 

double cell_properties::center_ny(unsigned int face_index) const {

    assert(face_index < 4);
    return face_center_normal_y[face_index]; 
}

// 15. Area of a face f 

double cell_properties::S_f(unsigned int face_index) const {

    assert(face_index < 4);
    return surface_area[face_index]; 
}


// 16. Size of the cell

double cell_properties::h() const {

	return h_size; 
}
