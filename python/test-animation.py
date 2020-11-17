import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np

def init_animation():
    global line
    line, = ax.plot(x, np.zeros_like(x))
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-1,1)

def animate(i):
    line.set_ydata(np.sin(2*np.pi*i / 50)*np.sin(x))
    return line,

fig = plt.figure()
ax = fig.add_subplot(111)
x = np.linspace(0, 2*np.pi, 200)

ani = matplotlib.animation.FuncAnimation(fig, animate, init_func=init_animation, frames=100)
ani.save('animation.gif', writer='imagemagick', fps=30)



##########################################
# import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation


# def update_line(num, data, line):
#     line.set_data(data[..., :num])
#     return line,

# # Fixing random state for reproducibility
# np.random.seed(19680801)

# # Set up formatting for the movie files
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

# fig1 = plt.figure()

# data = np.random.rand(2, 25)
# l, = plt.plot([], [], 'r-')
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.xlabel('x')
# plt.title('test')
# line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
#                                    interval=50, blit=True)
# line_ani.save('lines.mp4', writer=writer)

# fig2 = plt.figure()

# x = np.arange(-9, 10)
# y = np.arange(-9, 10).reshape(-1, 1)
# base = np.hypot(x, y)
# ims = []
# for add in np.arange(15):
#     ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))

# im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
#                                    blit=True)
# im_ani.save('im.mp4', writer=writer)
# ####################################################################
# # import numpy as np
# # import matplotlib.pyplot as plt
# # import matplotlib.animation as animation
# # from matplotlib import rcParams

# # # configure full path for ImageMagick
# # rcParams['animation.convert_path'] = r'/usr/bin/convert'

# # TWOPI = 2*np.pi

# # fig, ax = plt.subplots()

# # t = np.arange(0.0, TWOPI, 0.001)
# # s = np.sin(t)
# # l = plt.plot(t, s)

# # ax = plt.axis([0,TWOPI,-1,1])

# # redDot, = plt.plot([0], [np.sin(0)], 'ro')

# # def animate(i):
# #     redDot.set_data(i, np.sin(i))
# #     return redDot,

# # # create animation using the animate() function with no repeat
# # myAnimation = animation.FuncAnimation(fig, animate, frames=np.arange(0.0, TWOPI, 0.1), \
# #                                       interval=10, blit=True, repeat=True)
# # plt.show()
# # save animation at 30 frames per second

# # myAnimation.save('myAnimation.gif', writer='imagemagick', fps=30)
# # # myAnimation.save('myAnimation.gif', writer=animation.PillowWriter, fps=30)
# # writergif = animation.PillowWriter(fps=30)
# # myAnimation.save('testgif.gif', writer=writergif)

# ############################################################

# # import numpy as np
# # import matplotlib as mpl
# # import matplotlib.pyplot as plt
# # import matplotlib.animation as animation

# # mpl.use('Agg')
# # # mpl.rcParams['animation.ffmpeg_path'] = "C:/Path/To/Image/Magick/ffmpeg.exe"
# # # For Imagemagick 7.0, convert.exe is replaced by magick.exe
# # mpl.rcParams['animation.convert_path'] = r'/usr/bin/magick'

# # # Save the animation
# # ani = animation.FuncAnimation(fig, animateFunction, frames=N, interval=2, blit=True, repeat=False)
# # ani.save("./output.gif", writer='imagemagick', fps=60, bitrate=-1)
# # plt.show()

# #######################################################################

# # from matplotlib import animation
# # from matplotlib.animation import FuncAnimation, PillowWriter

# # # First set up the figure, the axis, and the plot element we want to animate
# # fig = plt.figure()
# # ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
# # line, = ax.plot([], [], lw=2)

# # # initialization function: plot the background of each frame
# # def init():
# #     line.set_data([], [])
# #     return line,

# # # animation function.  This is called sequentially
# # def animate(i):
# #     x = np.linspace(0, 2, 1000)
# #     y = np.sin(2 * np.pi * (x - 0.01 * i))
# #     line.set_data(x, y)
# #     return line,
# # # call the animator.  blit=True means only re-draw the parts that have changed.
# # anim = animation.FuncAnimation(fig, animate, init_func=init,frames=200, interval=20, blit=True)

# # save the animation as an mp4.  This requires ffmpeg or mencoder to be
# # installed.  The extra_args ensure that the x264 codec is used, so that
# # the video can be embedded in html5.  You may need to adjust this for
# # your system: for more information, see
# # http://matplotlib.sourceforge.net/api/animation_api.html
# # writergif = animation.PillowWriter(fps=30)
# # anim.save('basic_animation.gif', writer =PillowWriter(fps=25))# fps=30, extra_args=['-vcodec', 'libx264'])


# # fig, ax = plt.subplots()
# # x, ysin, ycos = [], [], []
# # ln1, = plt.plot([], [], 'ro')
# # ln2, = plt.plot([], [], 'm*')


# # def init():
# #     ax.set_xlim(0, 2*np.pi)
# #     ax.set_ylim(-1, 1)

# # def update(i):
# #     x.append(i)
# #     ysin.append(np.sin(i))
# #     ycos.append(np.cos(i))
# #     ln1.set_data(x, ysin)
# #     ln2.set_data(x, ycos)

# # ani = FuncAnimation(fig, update, np.linspace(0, 2*np.pi, 64), init_func=init)


# # writer = PillowWriter(fps=25)
# # ani.save("demo_sine.gif", writer=writer)