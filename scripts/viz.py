#!/usr/bin/env python3
import matplotlib.pyplot as plt
import os
import sys

# Add path handling to work with data directory
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
OUTPUT_DIR = os.path.join(DATA_DIR, 'output')

# Create directories if they don't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Данные
edge = [822, 441]
faces = [
    [441, 842, 822],
    [441, 822, 420]
]
vertices = {
    "441": [-0.070711, 0, -0.1],  # 441
    "842": [-0.06364, 0.007071, -0.1],  # 842
    "822": [-0.06364, 0.007071, -0.09],  # 822
    "420": [-0.070711, 0, -0.09]  # 420
}

# Функция для получения координат вершин по индексам
def get_vertex_coords(vertex_indices):
    return [vertices[str(i)] for i in vertex_indices]

# Получение координат вершин для каждого треугольника
face1_coords = get_vertex_coords(faces[0])
face2_coords = get_vertex_coords(faces[1])

# Преобразование координат в формат, пригодный для matplotlib
face1_x = [v[0] for v in face1_coords]
face1_y = [v[1] for v in face1_coords]
face1_z = [v[2] for v in face1_coords]

face2_x = [v[0] for v in face2_coords]
face2_y = [v[1] for v in face2_coords]
face2_z = [v[2] for v in face2_coords]

# Получение координат общего ребра
edge_coords = get_vertex_coords(edge)
edge_x = [v[0] for v in edge_coords]
edge_y = [v[1] for v in edge_coords]
edge_z = [v[2] for v in edge_coords]

# Создание фигуры и осей
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Рисование треугольников
ax.plot_trisurf(face1_x, face1_y, face1_z, triangles=[[0, 1, 2]], color='r', alpha=0.5)
ax.plot_trisurf(face2_x, face2_y, face2_z, triangles=[[0, 1, 2]], color='b', alpha=0.5)

# Рисование общего ребра
ax.plot(edge_x, edge_y, edge_z, color='g', linewidth=2)
handles = []
legend_labels = []
# Отображение всех точек и подпись их координат
for vertex, index in zip(vertices.values(), vertices.keys()):
    scatter = ax.scatter(vertex[0], vertex[1], vertex[2], color='k', s=50)
    ax.text(vertex[0], vertex[1], vertex[2]+0.001, f"({index})", fontsize=11, color='r')
    handles.append(scatter)
    legend_labels.append(f"{index}: {vertex[0]}, {vertex[1]}, {vertex[2]}")
     
# Добавление легенды
ax.legend(handles, legend_labels, title="Vertex Indices and Coordinates")

# Настройка осей
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Отображение графика
plt.show()