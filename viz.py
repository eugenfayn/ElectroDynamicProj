import matplotlib.pyplot as plt

def read_shared_edges(filename):
    edges = []
    faces = []
    vertices = []

    with open(filename, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith("Edge:"):
                edge = list(map(int, lines[i].split()[1:]))
                edges.append(edge)
                i += 1
                while i < len(lines) and not lines[i].startswith("Edge:"):
                    if lines[i].startswith("Face"):
                        face = list(map(int, lines[i].split()[2:]))
                        faces.append(face)
                        i += 1
                        vertex_coords = []
                        while i < len(lines) and not lines[i].startswith("Face") and not lines[i].startswith("Edge:"):
                            vertex = list(map(float, lines[i].split()))
                            vertex_coords.append(vertex)
                            i += 1
                        vertices.append(vertex_coords)
            else:
                i += 1
    return edges, faces, vertices

def plot_edges_and_faces(edges, faces, vertices):
    fig, ax = plt.subplots()
    for i in range(0):
        edge = edges[i]
        face1 = faces[2 * i]
        face2 = faces[2 * i + 1]
        vertex_coords1 = vertices[2 * i]
        vertex_coords2 = vertices[2 * i + 1]

        # Plot the shared edge
        ax.plot([vertex_coords1[0][0], vertex_coords1[1][0]], [vertex_coords1[0][1], vertex_coords1[1][1]], 'r-')

        # Plot the first face
        ax.plot([vertex_coords1[0][0], vertex_coords1[2][0]], [vertex_coords1[0][1], vertex_coords1[2][1]], 'b-')
        ax.plot([vertex_coords1[1][0], vertex_coords1[2][0]], [vertex_coords1[1][1], vertex_coords1[2][1]], 'b-')

        # Plot the second face
        ax.plot([vertex_coords2[0][0], vertex_coords2[2][0]], [vertex_coords2[0][1], vertex_coords2[2][1]], 'g-')
        ax.plot([vertex_coords2[1][0], vertex_coords2[2][0]], [vertex_coords2[1][1], vertex_coords2[2][1]], 'g-')

    ax.set_aspect('equal')
    plt.show()

edges, faces, vertices = read_shared_edges('shared_edges.txt')
plot_edges_and_faces(edges, faces, vertices)
