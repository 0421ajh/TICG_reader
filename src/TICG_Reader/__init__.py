#Reader for TICG simulation data in OVITO.
#
#+-----------------------------------------------+
#| By:                                           |
#|     - Jihun Ahn                               |
#| Date:                                         |
#|     - 2024.07.16                              |
#|                                               |
#+-----------------------------------------------+


from ovito.data import DataCollection
from ovito.io import FileReaderInterface, import_file
from typing import Callable, Any
import re

class TICG_Reader(FileReaderInterface):
    @staticmethod
    def detect(filename: str):
        try:
            with open(filename, "r") as f:
                _    = f.readline()
                line = f.readline()
                return line.strip() == "MC simulation of coarse grain block copolymer"
            #check their is PSF file
            if not os.path.exists('*.psf'):
                return False
        except OSError:
            return False


    def scan(self, filename: str, register_frame: Callable[..., None]):
        expr      = re.compile(r'^\d+$')
        label_num = 1
        with open(filename, "r") as f:
            for line_number, line in enumerate(f):
                match = re.match(expr, line.strip())
                if match:
                    num_particles = int(line.strip())
                    label = f"Frame {label_num}"
                    label_num += 1
                    register_frame(frame_info=(line_number+1, num_particles), label=label)
    def parse(self, data: DataCollection, filename: str, frame_info: tuple[int, int], **kwargs: Any):
        starting_line_number, num_particles = frame_info

        with open(filename, "r") as f:
            for _ in range(starting_line_number + 1):
                f.readline()

            particles = data.create_particles(count=num_particles)
            positions = particles.create_property("Position")

            for i in range(num_particles):
                parts = f.readline().strip().split()
                x, y, z = map(float, parts[1:4])
                positions[i] = [x, y, z]
                #  positions[i] = [float(coord) for coord in f.readline().strip().split()]


if __name__ == "__main__":
    pipeline = import_file("bead_out.xyz", input_format=TICG_Reader)
    for frame in range(pipeline.source.num_frames):
        data = pipeline.compute(frame)

