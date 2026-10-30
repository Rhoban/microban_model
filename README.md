# microban_model

The model of Microban v2 (the Radxa version), shared by [microban](https://github.com/Rhoban/microban) (the robot's
code), [microban_sim](https://github.com/Rhoban/microban_sim) (Rhoban's robot simulator, which `make sim` runs) and the
simulators built on it (the hackathon's [climbing simulator](https://github.com/rock-climbing-robots/simulator)). MuJoCo XML and URDF, in two versions of the same robot that differ
only by the forearms and the feet:

| | MuJoCo | URDF (kinematics and inertias) | mass |
|---|---|---|---|
| walking forearms and feet (default) | `mjcf/robot.xml` | `urdf/robot.urdf` | 769.8 g |
| climbing: ball hands, climbing feet | `mjcf/robot_climb.xml` | `urdf/robot_climb.urdf` | 799.3 g (cut feet, below) |

```python
import microban_model
microban_model.mjcf()               # path of mjcf/robot.xml; climb=True: robot_climb.xml
microban_model.urdf(climb=True)     # path of urdf/robot_climb.urdf
microban_model.JOINTS               # the 19 actuated joints, in the models' order
```

Install: `uv add "microban-model @ git+https://github.com/Rhoban/microban_model"` (the meshes are stored with Git LFS:
uv fetches them with `[tool.uv.sources] microban-model = { git = "...", lfs = true }`).

## This fork: the climbing feet cut

A fork of [Rhoban/microban_model](https://github.com/Rhoban/microban_model) at 8e42ef7, for the rock-climbing-robots
hackathon team whose robot has its climbing feet cut: the end of each foot's bottom plate is cut off (a slanted cut,
y 21 to 48 mm and z -13.5 to -1.4 mm in the foot mesh's frame). Only the climbing feet changed, on both feet, in
`robot_climb.xml` and `robot_climb.urdf`; the walking model, the hands and every other part are Rhoban's.

* `mjcf/assets/foot__configuration_right.stl`: the cut foot, in the same frame and units as before, so the model
  places it as it did.
* `mjcf/assets/foot__configuration_right_00.stl` ... `_50.stl`: its collision pieces, 51 convex pieces (CoACD,
  threshold 0.02, at most 64 vertices each) instead of 24; they follow the foot's surface to 0.6 mm.
* The feet's inertials (MJCF and URDF): Microban's CAD values minus the removed part (6.57 cm3 at the foot's mean
  density, 645 kg/m3): 29.98 g -> 25.74 g per foot, the centre of mass 2.7 mm closer to the ankle.

**To change the foot again** (a new cut, a new CAD export): on this branch, with git-lfs installed,

```bash
uv run --with numpy --with coacd==1.0.14 tools/cut_foot.py NEW_FOOT.stl            # in Rhoban's foot frame
uv run --with numpy --with coacd==1.0.14 tools/cut_foot.py NEW_FOOT.stl --rotate-x180  # as exported from Onshape
```

It writes the drawn foot, its collision pieces and both feet's inertials (from Rhoban's commit and the new mesh;
its docstring says how); run on the current foot it reproduces this branch exactly. Then commit, push, and pin the new
commit where the model is used (the team's repository: `pyproject.toml` and `MODEL_COMMIT`).

The STL files remain under the CAD licence (CC BY-NC-SA 4.0, Copyright (c) 2026 Marc Duclusaud, see `LICENSE`); the
changed foot and its pieces are shared under the same licence.

## Conventions

* **Every joint's zero is its motor's zero** (Dynamixel Wizard: 180 deg). The shoulder pitch at 0 holds the arm
  forward, horizontal; along the body it is +pi/2.
* Joint names are microban's (`left_knee`, ...); the trunk frame: x forward, y left, z up. The IMU is the site `imu`
  (sensors `orientation`, `angular-velocity`, `imu_accel`...); the feet frames are the sites `left_foot`, `right_foot`;
  the climbing version adds the 25 mm balls' centres `left_hand_ball`, `right_hand_ball`, `left_foot_ball`,
  `right_foot_ball`.

## What is in it

* Geometry and masses from the Onshape CAD (onshape-to-robot 1.8.3); the walking forearms and feet are those of
  Microban v1.
* Actuators: XL330 as torque motors (class `xl330`: armature, damping, friction loss), for a servo model such as
  [BAM](https://github.com/Rhoban/bam) that computes the torque.
* Collisions (group 3): boxes for the servos that stick out, the trunk as four convex hulls, the other parts as convex
  pieces, primitives for the walking forearms and feet, the hand / foot necks and spheres for the balls of the
  climbing version. Hands and feet collide with the rest of the robot (contype 3 / conaffinity 2), the other parts
  only with the world.
* No scene: each user adds its floor, lights or wall (MuJoCo's `MjSpec`, or an XML that includes the robot).
