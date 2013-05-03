c--------------------------------------------------------------------------
c   Tracelevel flags
c--------------------------------------------------------------------------

      integer instruction_trace
      parameter (instruction_trace = 1)
      integer contraction_trace
      parameter (contraction_trace = 2)
      integer proc_trace
      parameter (proc_trace = 4)

      common /trace/tracelevel, current_op, current_line,trace,
     *              call_marker, dryrun, simulator, pardo_timer,
     *              pardo_block_wait_timer,
     *              pardo_act_timer, pardo_tserver_timer,
     *              pardo_ovrhead_timer,
     *              pardo_times_timer,
     *              current_instr_timer,
     *              current_instr_blk_timer,
     *              current_instr_mpi_timer,
     *              current_instr_mpino_timer,
     *              current_instr_unit_timer,
     *              current_instr_allocate_timer,
     *              current_instr_total_timer
      integer tracelevel
      integer current_op
      integer current_line
      integer call_marker
      integer pardo_timer, pardo_block_wait_timer
      integer pardo_act_timer, pardo_tserver_timer, pardo_ovrhead_timer
      integer pardo_times_timer
      integer current_instr_timer, current_instr_blk_timer
      integer current_instr_mpi_timer
      integer current_instr_mpino_timer
      integer current_instr_unit_timer
      integer current_instr_allocate_timer
      integer current_instr_total_timer
      logical trace, dryrun, simulator
