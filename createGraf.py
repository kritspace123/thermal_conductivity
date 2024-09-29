import numpy as np
import plotly.graph_objs as go

# Задаем сетку
x = np.linspace(0, 1, 10)
t = np.linspace(0, 1, 100)
x, t = np.meshgrid(x, t)

# Определяем функцию z
z = -0.5 * (x**4) + (x**2) - x + t * x + 2 * (t**2) - t * np.exp(x)

# Создаем 3D поверхность
surface = go.Surface(z=z, x=x, y=t)

# Настраиваем макет с явными осями координат
layout = go.Layout(
    title='3D график с осями координат',
    scene=dict(
        xaxis=dict(
            title='X',
            showgrid=True,
            zeroline=True,
            showline=True,
            ticks='outside',
            showticklabels=True,
        ),
        yaxis=dict(
            title='t',
            showgrid=True,
            zeroline=True,
            showline=True,
            ticks='outside',
            showticklabels=True,
        ),
        zaxis=dict(
            title='Z',
            showgrid=True,
            zeroline=True,
            showline=True,
            ticks='outside',
            showticklabels=True,
        )
    )
)

# Создаем фигуру
fig = go.Figure(data=[surface], layout=layout)

# Визуализируем график
fig.show()
