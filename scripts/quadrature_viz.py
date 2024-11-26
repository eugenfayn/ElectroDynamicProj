import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os

def read_mesh_data(filename):
    data = {}
    
    with open(filename, 'r') as f:
        # Read vertices
        assert f.readline().strip() == "Vertices"
        n_vertices = int(f.readline())
        vertices = []
        for _ in range(n_vertices):
            x, y, z = map(float, f.readline().split())
            vertices.append([x, y, z])
        
        # Read faces
        assert f.readline().strip() == "Faces"
        n_faces = int(f.readline())
        faces = []
        for _ in range(n_faces):
            v1, v2, v3 = map(int, f.readline().split())
            faces.append([v1, v2, v3])
            
        # Read quadrature points
        assert f.readline().strip() == "QuadraturePoints"
        n_quad = int(f.readline())
        quad_points = []
        quad_weights = []
        for _ in range(n_quad):
            xi1, xi2, xi3, w = map(float, f.readline().split())
            quad_points.append([xi1, xi2, xi3])
            quad_weights.append(w)
    
    return {
        'vertices': np.array(vertices),
        'faces': np.array(faces),
        'quad_points': np.array(quad_points),
        'quad_weights': np.array(quad_weights)
    }

def barycentric_to_cartesian(bary_coords, triangle_vertices):
    """Convert barycentric coordinates to cartesian for a specific triangle"""
    result = np.zeros(3)
    for i in range(3):
        result += bary_coords[i] * triangle_vertices[i]
    return result

# Read the mesh data
data = read_mesh_data('mesh_data.txt')
vertices = data['vertices']
faces = data['faces']
quad_points = data['quad_points']
quad_weights = data['quad_weights']

# Print debug information
print("Number of vertices:", len(vertices))
print("Number of faces:", len(faces))
print("Number of quadrature points:", len(quad_points))
print("\nFirst few quadrature points:")
for i, (point, weight) in enumerate(zip(quad_points[:3], quad_weights[:3])):
    print(f"Point {i}: {point}, weight: {weight}")

# Create figure with two subplots
fig = plt.figure(figsize=(20, 10))

# First subplot: Full mesh
ax1 = fig.add_subplot(121, projection='3d')

# Plot full mesh
for face in faces:
    triangle_vertices = vertices[face]
    tri = Poly3DCollection([triangle_vertices])
    tri.set_alpha(0.3)
    tri.set_facecolor('blue')
    ax1.add_collection3d(tri)
    
    quad_points_cartesian = []
    for quad_point in quad_points:
        point = barycentric_to_cartesian(quad_point, triangle_vertices)
        quad_points_cartesian.append(point)
    
    quad_points_cartesian = np.array(quad_points_cartesian)
    ax1.scatter(quad_points_cartesian[:, 0], 
              quad_points_cartesian[:, 1], 
              quad_points_cartesian[:, 2], 
              color='red', s=50)

# Set limits for full mesh
x_min, x_max = vertices[:, 0].min(), vertices[:, 0].max()
y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()
z_min, z_max = vertices[:, 2].min(), vertices[:, 2].max()

margin = 0.1
x_range = x_max - x_min
y_range = y_max - y_min
z_range = z_max - z_min

ax1.set_xlim(x_min - margin * x_range, x_max + margin * x_range)
ax1.set_ylim(y_min - margin * y_range, y_max + margin * y_range)
ax1.set_zlim(z_min - margin * z_range, z_max + margin * z_range)

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.set_title('Full Mesh with Quadrature Points')

# Second subplot: First few triangles
ax2 = fig.add_subplot(122, projection='3d')

# Plot first few triangles
num_triangles_to_show = 3
for face_idx, face in enumerate(faces[:num_triangles_to_show]):
    triangle_vertices = vertices[face]
    tri = Poly3DCollection([triangle_vertices])
    tri.set_alpha(0.3)
    tri.set_facecolor(f'C{face_idx}')  # Different color for each triangle
    ax2.add_collection3d(tri)
    
    # Print vertices of this triangle
    print(f"\nTriangle {face_idx} vertices:")
    for i, vertex in enumerate(triangle_vertices):
        print(f"v{i+1}: {vertex}")
    
    quad_points_cartesian = []
    for quad_idx, quad_point in enumerate(quad_points):
        point = barycentric_to_cartesian(quad_point, triangle_vertices)
        quad_points_cartesian.append(point)
        print(f"Triangle {face_idx}, Quad point {quad_idx}: {point}")
    
    quad_points_cartesian = np.array(quad_points_cartesian)
    ax2.scatter(quad_points_cartesian[:, 0], 
               quad_points_cartesian[:, 1], 
               quad_points_cartesian[:, 2], 
               color=f'C{face_idx}', s=50, alpha=0.6)

# Set limits for detailed view
vertices_subset = vertices[np.unique(faces[:num_triangles_to_show])]
x_min, x_max = vertices_subset[:, 0].min(), vertices_subset[:, 0].max()
y_min, y_max = vertices_subset[:, 1].min(), vertices_subset[:, 1].max()
z_min, z_max = vertices_subset[:, 2].min(), vertices_subset[:, 2].max()

margin = 0.1
x_range = x_max - x_min
y_range = y_max - y_min
z_range = z_max - z_min

ax2.set_xlim(x_min - margin * x_range, x_max + margin * x_range)
ax2.set_ylim(y_min - margin * y_range, y_max + margin * y_range)
ax2.set_zlim(z_min - margin * z_range, z_max + margin * z_range)

ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.set_title(f'First {num_triangles_to_show} Triangles with Quadrature Points')

plt.tight_layout()
plt.show()
