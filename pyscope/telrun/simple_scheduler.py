from astroplan.scheduling import Scheduler, Scorer
from astroplan.utils import time_grid_from_range
from astroplan.constraints import AltitudeConstraint
from astropy import units as u

import numpy as np

class SimpleScheduler(Scheduler):
    """
    schedule blocks randomly
    """
    def __init__(self, *args, **kwargs):
        super(SimpleScheduler, self).__init__(*args, **kwargs)

    def _make_schedule(self, blocks):
        # gather all the constraints on each block into a single attribute
        for b in blocks:
            if b.constraints is None:
                b._all_constraints = self.constraints
            else:
                b._all_constraints = self.constraints + b.constraints

            # to make sure the Scorer has some constraint to work off of
            # and to prevent scheduling of targets below the horizon
            if b._all_constraints is None:
                b._all_constraints = [AltitudeConstraint(min=0*u.deg)]
                b.constraints = [AltitudeConstraint(min=0*u.deg)]
            elif not any(isinstance(c, AltitudeConstraint) for c in b._all_constraints):
                b._all_constraints.append(AltitudeConstraint(min=0*u.deg))
                if b.constraints is None:
                    b.constraints = [AltitudeConstraint(min=0*u.deg)]
                else:
                    b.constraints.append(AltitudeConstraint(min=0*u.deg))
            b.observer = self.observer

        # before we can schedule, we need to know where blocks meet the constraints
        scorer = Scorer(blocks, self.observer, self.schedule,
                        global_constraints=self.constraints)
        score_array = scorer.create_score_array(self.time_resolution)
        # now we have an array of the scores for the blocks at intervals of
        # ``time_resolution``. The scores range from zero to one, some blocks may have
        # higher scores than others, but we only care if they are greater than zero

        # we want to start from the beginning and start scheduling
        start_time = self.schedule.start_time
        current_time = start_time
        while current_time < self.schedule.end_time:
            scheduled = False
            i=0
            while i < len(blocks) and scheduled is False:
                block = blocks[i]
                # the schedule starts with only 1 slot
                if len(self.schedule.slots) == 1:
                    test_time = current_time
                # when a block is inserted, the number of slots increases
                else:
                    # a test transition between the last scheduled block and this one
                    transition = self.transitioner(schedule.observing_blocks[-1],
                                                   block, current_time, self.observer)
                    test_time = current_time + transition.duration
                # how many time intervals are we from the start
                start_idx = int((test_time - start_time)/self.time_resolution)
                duration_idx = int(block.duration/self.time_resolution)
                # if any score during the block's duration would be 0, reject it
                if any(score_array[i][start_idx:start_idx+duration_idx] == 0):
                    i +=1
                # if all of the scores are >0, accept and schedule it
                else:
                    if len(self.schedule.slots) >1:
                        self.schedule.insert_slot(current_time, transition)
                    self.schedule.insert_slot(test_time, block)
                    # advance the time and remove the block from the list
                    current_time = test_time + block.duration
                    scheduled = True
                    blocks.remove(block)
            # if every block failed, progress the time
            if i == len(blocks):
                current_time += self.gap_time
        return schedule