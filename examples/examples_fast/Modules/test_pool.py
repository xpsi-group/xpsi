""" Diagnose MPI parallel speed for the PT sampler.

Run from the Modules directory with the same launcher as the real job, e.g.

    mpiexec -n "$SLURM_NTASKS" python test_pool.py

It reports (1) how many MPI processes there are and on which nodes, (2) the
time of one likelihood call on every process, and (3) how long one batch of
likelihood calls takes through the MPI pool compared with the ideal time.
"""
import os
import sys
import time
import socket

# Optionally disable mpi4py's matched-probe receives, which hang with some
# Open MPI builds. Must be set before mpi4py.MPI is imported (also by xpsi).
#     MPROBE=0 mpiexec -n ... python test_pool.py
import mpi4py
if os.environ.get('MPROBE', '1') == '0':
    mpi4py.rc.recv_mprobe = False

import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def say(*args):
    if rank == 0:
        print(*args, flush=True)

say('MPI processes: %i' % size)
say('mpi4py %s, recv_mprobe=%s, OMPI_MCA_pml=%s'
    % (mpi4py.__version__, mpi4py.rc.recv_mprobe,
       os.environ.get('OMPI_MCA_pml')))
say('MPI library: %s' % MPI.Get_library_version().splitlines()[0])

# builds the model (and runs likelihood.check) on every process; the
# ``if __name__ == '__main__'`` block in main.py is not executed
t_setup = time.time()
from main import likelihood, prior
t_setup = time.time() - t_setup

import xpsi.Sample as S
from xpsi.Posterior import Posterior

S.posterior = Posterior(likelihood, prior)

# --- 1. where are the processes? -------------------------------------------
hosts = comm.gather(socket.gethostname(), root=0)
if rank == 0:
    from collections import Counter
    for h, n in sorted(Counter(hosts).items()):
        print('  %s: %i processes' % (h, n), flush=True)

# --- 2. serial timing of one likelihood call on every process ---------------
np.random.seed(1000 + rank)
pts = prior.draw(4)[0]
S.logl_func(pts[0]) # warm-up
t0 = time.time()
for q in pts[1:]:
    S.logp_func(q)
    S.logl_func(q)
t_call = (time.time() - t0) / (len(pts) - 1)

t_calls = comm.gather(t_call, root=0)
t_setups = comm.gather(t_setup, root=0)
if rank == 0:
    t_calls = np.array(t_calls)
    print('Model setup time per process: median %.1f s, max %.1f s'
          % (np.median(t_setups), np.max(t_setups)), flush=True)
    print('One likelihood call, serial on each process: '
          'min %.3f s, median %.3f s, max %.3f s'
          % (t_calls.min(), np.median(t_calls), t_calls.max()), flush=True)
    slow = np.where(t_calls > 3 * np.median(t_calls))[0]
    if len(slow):
        print('  slow ranks (>3x median):', slow.tolist(), flush=True)

comm.Barrier()

# --- 3. one batch through the pool, as ptemcee does it ----------------------
from schwimmbad import MPIPool

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    nworkers = pool.size
    nbatch = 2 * nworkers # two calls per worker
    batch = prior.draw(nbatch)[0]

    print('\nSending %i calls through the pool...' % nbatch, flush=True)
    done = [0]
    t0 = time.time()
    def progress(_):
        done[0] += 1
        if done[0] == 1 or done[0] % 50 == 0 or done[0] == nbatch:
            print('  %i / %i results back after %.2f s'
                  % (done[0], nbatch, time.time() - t0), flush=True)
    results = list(pool.map(S.logl_func, list(batch), callback=progress))
    t_pool = time.time() - t0

    ideal = nbatch * np.median(t_calls) / nworkers
    print('\nPool: %i calls on %i workers took %.2f s '
          '(ideal ~%.2f s, i.e. %.1fx slower than ideal)'
          % (nbatch, nworkers, t_pool, ideal, t_pool / ideal), flush=True)
    print('Finite results: %i / %i'
          % (np.sum(np.isfinite(results)), nbatch), flush=True)

    # what one ptemcee iteration would cost with the settings in main.py
    for nw, nt in [(256, 2)]:
        n = nw * nt
        print('Expected time per iteration for nwalkers=%i, ntemps=%i: '
              '~%.2f s' % (nw, nt, n * np.median(t_calls) / nworkers
                           * max(t_pool / ideal, 1.0)), flush=True)
