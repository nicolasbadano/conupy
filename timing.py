# timing.py encoding: utf-8
import time

class Timing(object):
    class Timer(object):
        def __init__(self, timing, name):
            self.timing = timing
            self.name = name

        def __enter__(self):
            self.start = time.clock()
            return self

        def __exit__(self, *args):
            ttime, num = self.timing.timers.get(self.name, (0, 0))
            self.timing.timers[self.name] = (ttime + time.clock() - self.start, num + 1)

    def __init__(self):
        self.timers = {}

    def getTimer(self, name):
        return Timing.Timer(self, name)

    def dump(self):
        for key, (ttime, num) in self.timers.iteritems():
            print "Code %s was called %i times and took an average of %f ms" % (key, num, ttime / num * 1000)
