import sys
import os
import customtkinter
import numpy as np
from mpi4py import MPI
from PIL import Image

def call_parallel_solver(refinement, timesteps, velocity, dt, viscosity):
    comm = MPI.COMM_SELF.Spawn(
        sys.executable,
        args = [os.path.dirname(os.path.abspath(__file__)) + '/solver.py', 
                refinement, timesteps, velocity, dt, viscosity],
        maxprocs=4
    )

    N = np.array(0, dtype='i')
    comm.Reduce(None, [N, MPI.INT], op=MPI.SUM, root=MPI.ROOT)

customtkinter.set_default_color_theme("dark-blue")

class OptionFrame(customtkinter.CTkFrame):
    def __init__(self, master, image_update):
        super().__init__(master)
        self.master = master
        self.image_update = image_update
        title_font = customtkinter.CTkFont(family="<family name>", weight='bold')

        #Define UI widgets in option menu
        self.title = customtkinter.CTkLabel(self, text='CFD Simulation', font=title_font)
        self.button = customtkinter.CTkButton(self, text="compute solution", command=self.button_callbck)
        self.timestepsL = customtkinter.CTkLabel(self, text="Timesteps:")
        self.timesteps = customtkinter.CTkEntry(self, placeholder_text="Enter timesteps")
        self.viscosityL = customtkinter.CTkLabel(self, text="Viscosity:")
        self.viscosity = customtkinter.CTkEntry(self, placeholder_text="Enter viscosity")
        self.refinementL = customtkinter.CTkLabel(self, text="Refinement:")
        self.refinement = customtkinter.CTkEntry(self, placeholder_text="Enter refinement")
        self.dtL = customtkinter.CTkLabel(self, text="Dt:")
        self.dt = customtkinter.CTkEntry(self, placeholder_text="Enter dt")
        self.velocityL = customtkinter.CTkLabel(self, text="Velocity:")
        self.velocity = customtkinter.CTkEntry(self, placeholder_text="Enter velocity")

        #Position widgets in grid
        self.title.grid(row=0, column=0, padx=20, pady=5, sticky="w")
        self.button.grid(row=1, column=0, padx=20, pady=5, sticky="w")
        self.timestepsL.grid(row=2, column=0, padx=20, pady=5, sticky="w")
        self.timesteps.grid(row=3, column=0, padx=20, pady=0, sticky="w")
        self.viscosityL.grid(row=4, column=0, padx=20, pady=5, sticky="w")
        self.viscosity.grid(row=5, column=0, padx=20, pady=0, sticky="w")
        self.refinementL.grid(row=6, column=0, padx=20, pady=5, sticky="w")
        self.refinement.grid(row=7, column=0, padx=20, pady=0, sticky="w")
        self.dtL.grid(row=8, column=0, padx=20, pady=5, sticky="w")
        self.dt.grid(row=9, column=0, padx=20, pady=0, sticky="w")
        self.velocityL.grid(row=10, column=0, padx=20, pady=5, sticky="w")
        self.velocity.grid(row=11, column=0, padx=20, pady=0, sticky="w")

    #Execute simulation given parameters
    def button_callbck(self):
        call_parallel_solver(self.refinement.get(), self.timesteps.get(), 
                           self.velocity.get(), self.dt.get(), self.viscosity.get())
        self.image_update()

class App(customtkinter.CTk):
    def __init__(self):
        #Initialise application UI
        super().__init__()
        customtkinter.set_appearance_mode("dark")
        self.title("CFD Solver")
        self.geometry("1000x600")
        self.image = customtkinter.CTkImage(dark_image=Image.open("stream_plot.png"), size=(640, 480))
        self.label = customtkinter.CTkLabel(master=self, image=self.image, text='')
        self.label.grid(row=0, column=1, padx=20, pady=20, sticky="w")
        self.option_frame = OptionFrame(self, self.image_update)
        self.option_frame.grid(row=0, column = 0, padx=20, pady=20, sticky="w")


    #Update image following simulation execution
    def image_update(self):
        self.image = customtkinter.CTkImage(dark_image=Image.open("stream_plot.png"), size=(640, 480))
        self.label.configure(image=self.image)

app = App()
app.mainloop()