<script lang="ts">
  import { onMount } from 'svelte';
  import { invoke } from '@tauri-apps/api/core';
  import { open } from '@tauri-apps/plugin-shell';
  import { load } from '@tauri-apps/plugin-store';
  import { currentView, pipelineStatus, logLines, reportPath, condaEnvPath } from '../stores';
  import FilePicker from '../components/FilePicker.svelte';
  import AdvancedOptions from '../components/AdvancedOptions.svelte';
  import type { PipelineConfig, ToolDef } from '../types';

  // ── Tool definitions ─────────────────────────────────────────────────────
  const TOOLS: ToolDef[] = [
    { id: 'run',            category: 'Simulation',  label: 'Run Pipeline',         description: 'Full end-to-end pipeline (prepare → simulate → sequence → map → metrics → report)' },
    { id: 'coverage-curve', category: 'Simulation',  label: 'Coverage Curve',       description: 'Sweep across CT / capture-fraction values and plot coverage depth curves' },
    { id: 'build-probes',   category: 'Probes',      label: 'Build Probes',         description: 'Design a probe set from target sequences' },
    { id: 'assess-probes',  category: 'Probes',      label: 'Assess Probes',        description: 'Combined probe coverage + cross-reactivity assessment' },
    { id: 'xreact',         category: 'Analysis',    label: 'Cross-Reactivity',     description: 'Check probe cross-reactivity against genomes or other probes' },
    { id: 'panel-qc',       category: 'Analysis',    label: 'Panel QC',             description: 'Assess target panel discriminability between species' },
    { id: 'identify',       category: 'Analysis',    label: 'Identify Species',     description: 'Call species present/absent/ambiguous from detection patterns' },
  ];

  const CATEGORIES = [...new Set(TOOLS.map(t => t.category))];

  let selectedTool = TOOLS[0].id;
  let runError = '';
  let launching = false;

  $: tool = TOOLS.find(t => t.id === selectedTool)!;

  // ── Shared state ──────────────────────────────────────────────────────────
  let outdir = '';
  let threads = '4';
  let outputPrefix = '';
  let reportMode = 'full';
  let cleanup = false;

  // ── run ───────────────────────────────────────────────────────────────────
  let r_targets = '';
  let r_probes = '';
  let r_distractors = '';
  let r_sample = '';
  let r_sampleIsFile = true;
  let r_sampleInline = '';
  let r_distractorMode: 'fraction' | 'ct' = 'fraction';
  let r_distractorFraction = '0.9';
  let r_ct = '25';
  let r_ctBaseline = '20';
  let r_ctBaselineFraction = '0.01';
  let r_simulateMode = 'thermodynamic';
  let r_numFragments = '10000';
  let r_captureFraction = '0.5';
  let r_hybTemp = '70';
  let r_readLength = '120';
  let r_seed = '';
  let r_genomes = '';
  let r_sampleTargetMap = '';
  let r_groups = '';
  let r_hostFasta = '';
  let r_identify = false;
  let r_runName = '';
  // fragment params
  let r_fragLenMean = '175';
  let r_fragLenMin = '150';
  let r_fragLenMax = '200';
  let r_minimapPreset = 'sr';
  // sequencing
  let r_numSequences = '';
  let r_outputFormat = 'fasta';
  let r_readSimulator = 'perfect';
  let r_sequencerProfile = '';
  let r_coverageDepth = '1.0';
  let r_pairedEnd = false;
  let r_peFargLenMean = '200';
  let r_peFargLenSd = '50';
  let r_longReadLenMean = '';
  let r_longReadLenSd = '';
  let r_badreadGlitches = '';
  let r_badreadJunkReads = '';
  let r_badreadRandomReads = '';
  let r_badreadChimeras = '';
  // CT calibration
  let r_ctEfficiency = '1.0';
  let r_ctCalibMode: 'single' | 'two-point' = 'single';
  let r_ctCal1 = '20,0.01';
  let r_ctCal2 = '25,0.003';
  // extra optional
  let r_distractorGroups = '';
  let r_identifyIdentThreshold = '90';
  let r_identifyMinUniq = '1';

  // ── build-probes ──────────────────────────────────────────────────────────
  let bp_targets = '';
  let bp_method = 'tile';
  let bp_probeLength = '120';
  let bp_step = '-60';
  let bp_catchStride = '60';
  let bp_catchMismatches = '5';
  let bp_catchExtension = '0';
  let bp_catchCoverage = '1.0';
  let bp_catchMinhashThreshold = '0.6';
  let bp_syottiMismatches = '40';
  let bp_syottiSeedLen = '20';
  let bp_ptStep = '1';
  let bp_ptIdentity = '0.9';
  let bp_ptCoverage = '0.9';
  let bp_ptBatchSize = '100';
  let bp_ptMaxPanelSize = '';
  let bp_ptMinDepth = '1';
  let bp_ptMaxIterations = '20';
  let bp_ptMinCoverageGain = '0.001';
  let bp_minGc = '0.20';
  let bp_maxGc = '0.80';
  let bp_maxNFrac = '0.05';
  let bp_noNInProbes = false;
  let bp_dustThreshold = '2.0';
  let bp_dustWindow = '64';
  let bp_maxMaskedFrac = '0.25';
  let bp_collapseThreshold = '0.95';
  let bp_dedupThreshold = '0.95';
  let bp_proximity = '50';
  let bp_minimapPreset = 'sr';
  let bp_threshold = '80';
  let bp_refineThreshold = '80';
  let bp_refineMode: 'none' | 'iterations' | 'stable' = 'none';
  let bp_refineIterations = '3';
  let bp_skipAssess = false;
  let bp_genomes = '';

  // ── assess-probes ─────────────────────────────────────────────────────────
  let ap_targets = '';
  let ap_probes = '';
  let ap_genomes = '';
  let ap_threshold = '80';
  let ap_noIndividual = false;
  let ap_proximity = '50';
  let ap_minimapPreset = 'sr';
  let ap_gapMinLength = '';
  let ap_refineThreshold = '80';
  let ap_refineMode: 'none' | 'iterations' | 'stable' = 'none';
  let ap_refineIterations = '3';

  // ── xreact ───────────────────────────────────────────────────────────────
  let xr_probes = '';
  let xr_against = '';
  let xr_self = false;
  let xr_threshold = '80';
  let xr_minimapPreset = 'sr';

  // ── panel-qc ─────────────────────────────────────────────────────────────
  let pq_targets = '';
  let pq_sampleTargetMap = '';
  let pq_identityThreshold = '90';
  let pq_minimapPreset = 'sr';

  // ── identify ─────────────────────────────────────────────────────────────
  let id_detectedDetail = '';
  let id_sampleTargetMap = '';
  let id_targetSimilarity = '';
  let id_targets = '';
  let id_identityThreshold = '90';
  let id_minUniqueTargets = '1';
  let id_minimapPreset = 'sr';

  // ── coverage-curve ────────────────────────────────────────────────────────
  let cc_targets = '';
  let cc_probes = '';
  let cc_distractors = '';
  let cc_sampleIsFile = true;
  let cc_sample = '';
  let cc_sampleInline = '';
  let cc_genomes = '';
  let cc_sampleTargetMap = '';
  // Simulation
  let cc_numFragments = '10000';
  let cc_simulateMode = 'thermodynamic';
  // Distractor mode: CT or fraction
  let cc_distractorMode: 'ct' | 'fraction' = 'ct';
  // CT sweep
  let cc_ctSweep = false;
  let cc_ctFixed = '25';
  let cc_ctList = '20 25 30';
  // CT calibration
  let cc_ctCalibMode: 'single' | 'two-point' = 'single';
  let cc_ctBaseline = '20';
  let cc_ctBaselineFraction = '0.01';
  let cc_ctEfficiency = '1.0';
  let cc_ctCal1 = '20,0.01';
  let cc_ctCal2 = '25,0.003';
  // Distractor fraction sweep (alternative to CT)
  let cc_distractorFracSweep = false;
  let cc_distractorFracFixed = '0.9';
  let cc_distractorFracList = '0.5 0.7 0.9';
  // Temperature sweep
  let cc_tempSweep = false;
  let cc_tempFixed = '70';
  let cc_tempList = '60 65 70 75';
  // Capture fraction sweep
  let cc_cfSweep = false;
  let cc_cfFixed = '0.5';
  let cc_cfList = '0.3 0.5 0.8';
  // Num sequences sweep
  let cc_nsSweep = false;
  let cc_nsFixed = '';
  let cc_nsList = '500 1000 5000';
  // Advanced
  let cc_seed = '';
  let cc_fragLenMean = '175';
  let cc_fragLenMin = '150';
  let cc_fragLenMax = '200';
  let cc_hostFasta = '';
  let cc_minimapPreset = 'sr';
  let cc_hostMinimapPreset = 'sr';

  // ── Persistence ──────────────────────────────────────────────────────────

  let _formStore: Awaited<ReturnType<typeof load>> | null = null;
  async function getFormStore() {
    if (!_formStore) _formStore = await load('form-state.json', { autoSave: false });
    return _formStore;
  }

  function captureState(toolId: string): Record<string, unknown> {
    const shared = { outdir, threads, outputPrefix, reportMode, cleanup };
    switch (toolId) {
      case 'run': return { ...shared,
        r_targets, r_probes, r_distractors, r_sampleIsFile, r_sample, r_sampleInline,
        r_distractorMode, r_distractorFraction, r_ct, r_ctBaseline, r_ctBaselineFraction,
        r_ctEfficiency, r_ctCalibMode, r_ctCal1, r_ctCal2,
        r_simulateMode, r_numFragments, r_captureFraction, r_hybTemp, r_readLength, r_seed,
        r_genomes, r_sampleTargetMap, r_groups, r_distractorGroups, r_hostFasta,
        r_identify, r_identifyIdentThreshold, r_identifyMinUniq, r_runName,
        r_fragLenMean, r_fragLenMin, r_fragLenMax, r_minimapPreset,
        r_numSequences, r_outputFormat, r_readSimulator, r_sequencerProfile, r_coverageDepth,
        r_pairedEnd, r_peFargLenMean, r_peFargLenSd,
        r_longReadLenMean, r_longReadLenSd,
        r_badreadGlitches, r_badreadJunkReads, r_badreadRandomReads, r_badreadChimeras,
      };
      case 'build-probes': return { ...shared,
        bp_targets, bp_method, bp_probeLength, bp_step,
        bp_catchStride, bp_catchMismatches, bp_catchExtension, bp_catchCoverage, bp_catchMinhashThreshold,
        bp_syottiMismatches, bp_syottiSeedLen,
        bp_ptStep, bp_ptIdentity, bp_ptCoverage, bp_ptBatchSize, bp_ptMaxPanelSize,
        bp_ptMinDepth, bp_ptMaxIterations, bp_ptMinCoverageGain,
        bp_minGc, bp_maxGc, bp_maxNFrac, bp_noNInProbes,
        bp_dustThreshold, bp_dustWindow, bp_maxMaskedFrac,
        bp_collapseThreshold, bp_dedupThreshold,
        bp_proximity, bp_minimapPreset, bp_threshold,
        bp_refineThreshold, bp_refineMode, bp_refineIterations,
        bp_skipAssess, bp_genomes,
      };
      case 'assess-probes': return { ...shared,
        ap_targets, ap_probes, ap_genomes, ap_threshold, ap_noIndividual,
        ap_proximity, ap_minimapPreset, ap_gapMinLength,
        ap_refineThreshold, ap_refineMode, ap_refineIterations,
      };
      case 'xreact': return { ...shared,
        xr_probes, xr_against, xr_self, xr_threshold, xr_minimapPreset,
      };
      case 'panel-qc': return { ...shared,
        pq_targets, pq_sampleTargetMap, pq_identityThreshold, pq_minimapPreset,
      };
      case 'identify': return { ...shared,
        id_detectedDetail, id_sampleTargetMap, id_targetSimilarity, id_targets,
        id_identityThreshold, id_minUniqueTargets, id_minimapPreset,
      };
      case 'coverage-curve': return { ...shared,
        cc_targets, cc_probes, cc_distractors, cc_sampleIsFile, cc_sample, cc_sampleInline,
        cc_genomes, cc_sampleTargetMap, cc_numFragments, cc_simulateMode,
        cc_distractorMode, cc_ctSweep, cc_ctFixed, cc_ctList,
        cc_ctCalibMode, cc_ctBaseline, cc_ctBaselineFraction, cc_ctEfficiency, cc_ctCal1, cc_ctCal2,
        cc_distractorFracSweep, cc_distractorFracFixed, cc_distractorFracList,
        cc_tempSweep, cc_tempFixed, cc_tempList,
        cc_cfSweep, cc_cfFixed, cc_cfList,
        cc_nsSweep, cc_nsFixed, cc_nsList,
        cc_seed, cc_fragLenMean, cc_fragLenMin, cc_fragLenMax,
        cc_hostFasta, cc_minimapPreset, cc_hostMinimapPreset,
      };
      default: return shared;
    }
  }

  function applyState(toolId: string, s: Record<string, unknown>) {
    const str = (k: string, d = '') => k in s && s[k] != null ? String(s[k]) : d;
    const bl  = (k: string, d = false) => k in s && s[k] != null ? Boolean(s[k]) : d;
    outdir = str('outdir'); threads = str('threads', '4');
    outputPrefix = str('outputPrefix'); reportMode = str('reportMode', 'full');
    cleanup = bl('cleanup');
    switch (toolId) {
      case 'run':
        r_targets = str('r_targets'); r_probes = str('r_probes'); r_distractors = str('r_distractors');
        r_sampleIsFile = bl('r_sampleIsFile', true); r_sample = str('r_sample'); r_sampleInline = str('r_sampleInline');
        r_distractorMode = str('r_distractorMode', 'fraction') as 'fraction' | 'ct';
        r_distractorFraction = str('r_distractorFraction', '0.9');
        r_ct = str('r_ct', '25'); r_ctBaseline = str('r_ctBaseline', '20');
        r_ctBaselineFraction = str('r_ctBaselineFraction', '0.01');
        r_ctEfficiency = str('r_ctEfficiency', '1.0');
        r_ctCalibMode = str('r_ctCalibMode', 'single') as 'single' | 'two-point';
        r_ctCal1 = str('r_ctCal1', '20,0.01'); r_ctCal2 = str('r_ctCal2', '25,0.003');
        r_simulateMode = str('r_simulateMode', 'thermodynamic');
        r_numFragments = str('r_numFragments', '10000'); r_captureFraction = str('r_captureFraction', '0.5');
        r_hybTemp = str('r_hybTemp', '70'); r_readLength = str('r_readLength', '120');
        r_seed = str('r_seed');
        r_genomes = str('r_genomes'); r_sampleTargetMap = str('r_sampleTargetMap');
        r_groups = str('r_groups'); r_distractorGroups = str('r_distractorGroups');
        r_hostFasta = str('r_hostFasta');
        r_identify = bl('r_identify');
        r_identifyIdentThreshold = str('r_identifyIdentThreshold', '90');
        r_identifyMinUniq = str('r_identifyMinUniq', '1');
        r_runName = str('r_runName');
        r_fragLenMean = str('r_fragLenMean', '175'); r_fragLenMin = str('r_fragLenMin', '150');
        r_fragLenMax = str('r_fragLenMax', '200'); r_minimapPreset = str('r_minimapPreset', 'sr');
        r_numSequences = str('r_numSequences'); r_outputFormat = str('r_outputFormat', 'fasta');
        r_readSimulator = str('r_readSimulator', 'perfect');
        r_sequencerProfile = str('r_sequencerProfile'); r_coverageDepth = str('r_coverageDepth', '1.0');
        r_pairedEnd = bl('r_pairedEnd');
        r_peFargLenMean = str('r_peFargLenMean', '200'); r_peFargLenSd = str('r_peFargLenSd', '50');
        r_longReadLenMean = str('r_longReadLenMean'); r_longReadLenSd = str('r_longReadLenSd');
        r_badreadGlitches = str('r_badreadGlitches'); r_badreadJunkReads = str('r_badreadJunkReads');
        r_badreadRandomReads = str('r_badreadRandomReads'); r_badreadChimeras = str('r_badreadChimeras');
        break;
      case 'build-probes':
        bp_targets = str('bp_targets'); bp_method = str('bp_method', 'tile');
        bp_probeLength = str('bp_probeLength', '120'); bp_step = str('bp_step', '-60');
        bp_catchStride = str('bp_catchStride', '60'); bp_catchMismatches = str('bp_catchMismatches', '5');
        bp_catchExtension = str('bp_catchExtension', '0'); bp_catchCoverage = str('bp_catchCoverage', '1.0');
        bp_catchMinhashThreshold = str('bp_catchMinhashThreshold', '0.6');
        bp_syottiMismatches = str('bp_syottiMismatches', '40'); bp_syottiSeedLen = str('bp_syottiSeedLen', '20');
        bp_ptStep = str('bp_ptStep', '1'); bp_ptIdentity = str('bp_ptIdentity', '0.9');
        bp_ptCoverage = str('bp_ptCoverage', '0.9'); bp_ptBatchSize = str('bp_ptBatchSize', '100');
        bp_ptMaxPanelSize = str('bp_ptMaxPanelSize');
        bp_ptMinDepth = str('bp_ptMinDepth', '1'); bp_ptMaxIterations = str('bp_ptMaxIterations', '20');
        bp_ptMinCoverageGain = str('bp_ptMinCoverageGain', '0.001');
        bp_minGc = str('bp_minGc', '0.20'); bp_maxGc = str('bp_maxGc', '0.80');
        bp_maxNFrac = str('bp_maxNFrac', '0.05'); bp_noNInProbes = bl('bp_noNInProbes');
        bp_dustThreshold = str('bp_dustThreshold', '2.0'); bp_dustWindow = str('bp_dustWindow', '64');
        bp_maxMaskedFrac = str('bp_maxMaskedFrac', '0.25');
        bp_collapseThreshold = str('bp_collapseThreshold', '0.95');
        bp_dedupThreshold = str('bp_dedupThreshold', '0.95');
        bp_proximity = str('bp_proximity', '50'); bp_minimapPreset = str('bp_minimapPreset', 'sr');
        bp_threshold = str('bp_threshold', '80');
        bp_refineThreshold = str('bp_refineThreshold', '80');
        bp_refineMode = str('bp_refineMode', 'none') as 'none' | 'iterations' | 'stable';
        bp_refineIterations = str('bp_refineIterations', '3');
        bp_skipAssess = bl('bp_skipAssess'); bp_genomes = str('bp_genomes');
        break;
      case 'assess-probes':
        ap_targets = str('ap_targets'); ap_probes = str('ap_probes'); ap_genomes = str('ap_genomes');
        ap_threshold = str('ap_threshold', '80'); ap_noIndividual = bl('ap_noIndividual');
        ap_proximity = str('ap_proximity', '50'); ap_minimapPreset = str('ap_minimapPreset', 'sr');
        ap_gapMinLength = str('ap_gapMinLength');
        ap_refineThreshold = str('ap_refineThreshold', '80');
        ap_refineMode = str('ap_refineMode', 'none') as 'none' | 'iterations' | 'stable';
        ap_refineIterations = str('ap_refineIterations', '3');
        break;
      case 'xreact':
        xr_probes = str('xr_probes'); xr_against = str('xr_against');
        xr_self = bl('xr_self'); xr_threshold = str('xr_threshold', '80');
        xr_minimapPreset = str('xr_minimapPreset', 'sr');
        break;
      case 'panel-qc':
        pq_targets = str('pq_targets'); pq_sampleTargetMap = str('pq_sampleTargetMap');
        pq_identityThreshold = str('pq_identityThreshold', '90');
        pq_minimapPreset = str('pq_minimapPreset', 'sr');
        break;
      case 'identify':
        id_detectedDetail = str('id_detectedDetail'); id_sampleTargetMap = str('id_sampleTargetMap');
        id_targetSimilarity = str('id_targetSimilarity'); id_targets = str('id_targets');
        id_identityThreshold = str('id_identityThreshold', '90');
        id_minUniqueTargets = str('id_minUniqueTargets', '1');
        id_minimapPreset = str('id_minimapPreset', 'sr');
        break;
      case 'coverage-curve':
        cc_targets = str('cc_targets'); cc_probes = str('cc_probes'); cc_distractors = str('cc_distractors');
        cc_sampleIsFile = bl('cc_sampleIsFile', true); cc_sample = str('cc_sample'); cc_sampleInline = str('cc_sampleInline');
        cc_genomes = str('cc_genomes'); cc_sampleTargetMap = str('cc_sampleTargetMap');
        cc_numFragments = str('cc_numFragments', '10000');
        cc_simulateMode = str('cc_simulateMode', 'thermodynamic');
        cc_distractorMode = str('cc_distractorMode', 'ct') as 'ct' | 'fraction';
        cc_ctSweep = bl('cc_ctSweep'); cc_ctFixed = str('cc_ctFixed', '25'); cc_ctList = str('cc_ctList', '20 25 30');
        cc_ctCalibMode = str('cc_ctCalibMode', 'single') as 'single' | 'two-point';
        cc_ctBaseline = str('cc_ctBaseline', '20'); cc_ctBaselineFraction = str('cc_ctBaselineFraction', '0.01');
        cc_ctEfficiency = str('cc_ctEfficiency', '1.0');
        cc_ctCal1 = str('cc_ctCal1', '20,0.01'); cc_ctCal2 = str('cc_ctCal2', '25,0.003');
        cc_distractorFracSweep = bl('cc_distractorFracSweep');
        cc_distractorFracFixed = str('cc_distractorFracFixed', '0.9');
        cc_distractorFracList = str('cc_distractorFracList', '0.5 0.7 0.9');
        cc_tempSweep = bl('cc_tempSweep'); cc_tempFixed = str('cc_tempFixed', '70'); cc_tempList = str('cc_tempList', '60 65 70 75');
        cc_cfSweep = bl('cc_cfSweep'); cc_cfFixed = str('cc_cfFixed', '0.5'); cc_cfList = str('cc_cfList', '0.3 0.5 0.8');
        cc_nsSweep = bl('cc_nsSweep'); cc_nsFixed = str('cc_nsFixed'); cc_nsList = str('cc_nsList', '500 1000 5000');
        cc_seed = str('cc_seed');
        cc_fragLenMean = str('cc_fragLenMean', '175'); cc_fragLenMin = str('cc_fragLenMin', '150');
        cc_fragLenMax = str('cc_fragLenMax', '200');
        cc_hostFasta = str('cc_hostFasta'); cc_minimapPreset = str('cc_minimapPreset', 'sr');
        cc_hostMinimapPreset = str('cc_hostMinimapPreset', 'sr');
        break;
    }
  }

  async function switchTool(toolId: string) {
    try {
      const store = await getFormStore();
      await store.set(`formState_${selectedTool}`, captureState(selectedTool));
      await store.set('lastTool', toolId);
      await store.save();
      const saved = await store.get<Record<string, unknown>>(`formState_${toolId}`);
      if (saved) applyState(toolId, saved);
    } catch { /* ignore store errors */ }
    selectedTool = toolId;
    runError = '';
  }

  onMount(async () => {
    try {
      const store = await getFormStore();
      const lastTool = await store.get<string>('lastTool');
      if (lastTool && TOOLS.some(t => t.id === lastTool)) selectedTool = lastTool;
      const saved = await store.get<Record<string, unknown>>(`formState_${selectedTool}`);
      if (saved) applyState(selectedTool, saved);
    } catch { /* first run or store error — use defaults */ }
  });

  // ── Build config and run ──────────────────────────────────────────────────
  function buildArgs(): Record<string, string> | null {
    const a: Record<string, string> = {};
    const add = (flag: string, val: string | number) => { const s = String(val); if (s !== '' && s !== 'undefined' && s !== 'null') a[flag] = s; };
    const flag = (f: string, on: boolean) => { if (on) a[f] = ''; };

    // Shared — outdir + output-prefix apply to all tools
    if (outdir) add('--outdir', outdir);
    if (outputPrefix) add('--output-prefix', outputPrefix);
    // report and cleanup are not supported by 'identify'
    if (selectedTool !== 'identify') {
      if (reportMode !== 'full') add('--report', reportMode);
      flag('--cleanup', cleanup);
    }

    switch (selectedTool) {
      case 'run': {
        if (!r_targets || !r_probes || !r_distractors) {
          runError = 'Targets, probes, and at least one distractor file are required.';
          return null;
        }
        add('--targets', r_targets);
        add('--probes', r_probes);
        add('--distractors', r_distractors);
        // sample
        const sampleVal = r_sampleIsFile ? r_sample : r_sampleInline;
        if (sampleVal) add('--sample', sampleVal);
        // distractor mode
        if (r_distractorMode === 'fraction') {
          add('--distractor-fraction', r_distractorFraction);
        } else {
          add('--ct', r_ct);
          if (r_ctCalibMode === 'two-point') {
            a['--ct-calibration'] = `${r_ctCal1}\t${r_ctCal2}`;
          } else {
            if (r_ctBaseline !== '20') add('--ct-baseline', r_ctBaseline);
            if (r_ctBaselineFraction !== '0.01') add('--ct-baseline-fraction', r_ctBaselineFraction);
            if (r_ctEfficiency !== '1.0') add('--ct-efficiency', r_ctEfficiency);
          }
        }
        add('--simulate-mode', r_simulateMode);
        add('--num-fragments', r_numFragments);
        add('--capture-fraction', r_captureFraction);
        if (r_simulateMode === 'thermodynamic') add('--hybridization-temperature', r_hybTemp);
        add('--read-length', r_readLength);
        if (r_seed) add('--seed', r_seed);
        add('--threads', threads);
        // sequencing
        if (r_numSequences) add('--num-sequences', r_numSequences);
        if (r_outputFormat !== 'fasta') add('--output-format', r_outputFormat);
        if (r_readSimulator !== 'perfect') {
          add('--read-simulator', r_readSimulator);
          if (r_sequencerProfile) add('--sequencer-profile', r_sequencerProfile);
          if (r_coverageDepth !== '1.0') add('--coverage-depth', r_coverageDepth);
          if (r_readSimulator === 'art') {
            flag('--paired-end', r_pairedEnd);
            if (r_pairedEnd) {
              add('--pe-frag-len-mean', r_peFargLenMean);
              add('--pe-frag-len-sd', r_peFargLenSd);
            }
          }
          if (r_readSimulator === 'badread') {
            if (r_longReadLenMean) add('--long-read-length-mean', r_longReadLenMean);
            if (r_longReadLenSd) add('--long-read-length-sd', r_longReadLenSd);
            if (r_badreadGlitches) add('--badread-glitches', r_badreadGlitches);
            if (r_badreadJunkReads) add('--badread-junk-reads', r_badreadJunkReads);
            if (r_badreadRandomReads) add('--badread-random-reads', r_badreadRandomReads);
            if (r_badreadChimeras) add('--badread-chimeras', r_badreadChimeras);
          }
        }
        // optional inputs
        if (r_genomes) add('--genomes', r_genomes);
        if (r_sampleTargetMap) add('--sample-target-map', r_sampleTargetMap);
        if (r_groups) add('--groups', r_groups);
        if (r_distractorGroups) add('--distractor-groups', r_distractorGroups);
        if (r_hostFasta) add('--host-fasta', r_hostFasta);
        flag('--identify', r_identify);
        if (r_identify) {
          if (r_identifyIdentThreshold !== '90') add('--identity-threshold', r_identifyIdentThreshold);
          if (r_identifyMinUniq !== '1') add('--min-unique-targets', r_identifyMinUniq);
        }
        if (r_runName) add('--run-name', r_runName);
        // fragment params
        add('--fragment-length-mean', r_fragLenMean);
        add('--fragment-length-min', r_fragLenMin);
        add('--fragment-length-max', r_fragLenMax);
        if (r_minimapPreset) add('--minimap-preset', r_minimapPreset);
        break;
      }

      case 'build-probes': {
        if (!bp_targets) { runError = 'Targets file is required.'; return null; }
        add('--targets', bp_targets);
        add('--method', bp_method);
        add('--probe-length', bp_probeLength);
        if (bp_method === 'tile') add('--step', bp_step);
        if (bp_method === 'catch-lite' || bp_method === 'catch') {
          add('--catch-probe-stride', bp_catchStride);
          add('--catch-mismatches', bp_catchMismatches);
          if (bp_catchExtension !== '0') add('--catch-extension', bp_catchExtension);
          if (bp_catchCoverage !== '1.0') add('--catch-coverage', bp_catchCoverage);
          if (bp_catchMinhashThreshold !== '0.6') add('--catch-minhash-threshold', bp_catchMinhashThreshold);
        }
        if (bp_method === 'syotti-lite') {
          add('--syotti-mismatches', bp_syottiMismatches);
          add('--syotti-seed-len', bp_syottiSeedLen);
        }
        if (bp_method === 'probetools-lite') {
          if (bp_ptStep !== '1') add('--pt-step', bp_ptStep);
          if (bp_ptIdentity !== '0.9') add('--pt-identity', bp_ptIdentity);
          if (bp_ptCoverage !== '0.9') add('--pt-coverage', bp_ptCoverage);
          if (bp_ptBatchSize !== '100') add('--pt-batch-size', bp_ptBatchSize);
          if (bp_ptMaxPanelSize) add('--pt-max-panel-size', bp_ptMaxPanelSize);
          if (bp_ptMinDepth !== '1') add('--pt-min-depth', bp_ptMinDepth);
          if (bp_ptMaxIterations !== '20') add('--pt-max-iterations', bp_ptMaxIterations);
          if (bp_ptMinCoverageGain !== '0.001') add('--pt-min-coverage-gain', bp_ptMinCoverageGain);
        }
        add('--min-gc', bp_minGc);
        add('--max-gc', bp_maxGc);
        if (bp_maxNFrac !== '0.05') add('--max-n-frac', bp_maxNFrac);
        flag('--no-n-in-probes', bp_noNInProbes);
        add('--dust-threshold', bp_dustThreshold);
        if (bp_dustWindow !== '64') add('--dust-window', bp_dustWindow);
        if (bp_maxMaskedFrac !== '0.25') add('--max-masked-frac', bp_maxMaskedFrac);
        if (bp_collapseThreshold !== '0.95') add('--collapse-threshold', bp_collapseThreshold);
        if (bp_dedupThreshold !== '0.95') add('--dedup-threshold', bp_dedupThreshold);
        if (bp_proximity !== '50') add('--proximity', bp_proximity);
        if (bp_minimapPreset !== 'sr') add('--minimap-preset', bp_minimapPreset);
        if (bp_threshold !== '80') add('--threshold', bp_threshold);
        if (bp_refineMode !== 'none') {
          add('--refine-threshold', bp_refineThreshold);
          if (bp_refineMode === 'iterations') add('--refine-iterations', bp_refineIterations);
          else flag('--refine-until-stable', true);
        }
        if (bp_genomes) add('--genomes', bp_genomes);
        flag('--skip-assess', bp_skipAssess);
        add('--threads', threads);
        break;
      }

      case 'assess-probes': {
        if (!ap_targets || !ap_probes) { runError = 'Targets and probes files are required.'; return null; }
        add('--targets', ap_targets);
        add('--probes', ap_probes);
        if (ap_genomes) add('--genomes', ap_genomes);
        if (ap_threshold !== '80') add('--threshold', ap_threshold);
        flag('--no-individual-targets', ap_noIndividual);
        if (ap_proximity !== '50') add('--proximity', ap_proximity);
        if (ap_minimapPreset !== 'sr') add('--minimap-preset', ap_minimapPreset);
        if (ap_gapMinLength) add('--gap-min-length', ap_gapMinLength);
        if (ap_refineMode !== 'none') {
          add('--refine-threshold', ap_refineThreshold);
          if (ap_refineMode === 'iterations') add('--refine-iterations', ap_refineIterations);
          else flag('--refine-until-stable', true);
        }
        add('--threads', threads);
        break;
      }

      case 'xreact': {
        if (!xr_probes) { runError = 'Probes file is required.'; return null; }
        if (!xr_against && !xr_self) {
          runError = 'Specify at least one of: genome FASTA(s) or Self cross-reactivity check.';
          return null;
        }
        add('--probes', xr_probes);
        if (xr_against) add('--against', xr_against);
        flag('--self', xr_self);
        if (xr_threshold !== '80') add('--threshold', xr_threshold);
        if (xr_minimapPreset !== 'sr') add('--minimap-preset', xr_minimapPreset);
        break;
      }

      case 'panel-qc': {
        if (!pq_targets || !pq_sampleTargetMap) { runError = 'Targets and sample-target-map are required.'; return null; }
        add('--targets', pq_targets);
        add('--sample-target-map', pq_sampleTargetMap);
        if (pq_identityThreshold !== '90') add('--identity-threshold', pq_identityThreshold);
        if (pq_minimapPreset !== 'sr') add('--minimap-preset', pq_minimapPreset);
        break;
      }

      case 'identify': {
        if (!id_detectedDetail || !id_sampleTargetMap) { runError = 'detected_detail.tsv and sample-target-map are required.'; return null; }
        add('--detected-detail', id_detectedDetail);
        add('--sample-target-map', id_sampleTargetMap);
        if (id_targetSimilarity) add('--target-similarity', id_targetSimilarity);
        else if (id_targets) add('--targets', id_targets);
        if (id_identityThreshold !== '90') add('--identity-threshold', id_identityThreshold);
        if (id_minUniqueTargets !== '1') add('--min-unique-targets', id_minUniqueTargets);
        if (id_minimapPreset !== 'sr') add('--minimap-preset', id_minimapPreset);
        break;
      }

      case 'coverage-curve': {
        if (!cc_targets || !cc_probes || !cc_distractors) { runError = 'Targets, probes, and distractors are required.'; return null; }
        add('--targets', cc_targets);
        add('--probes', cc_probes);
        add('--distractors', cc_distractors);
        if (cc_genomes) add('--genomes', cc_genomes);
        if (cc_sampleTargetMap) add('--sample-target-map', cc_sampleTargetMap);
        // Sample
        if (cc_sampleIsFile) {
          if (cc_sample) add('--sample', cc_sample);
        } else if (cc_sampleInline.trim()) {
          a['--sample'] = cc_sampleInline.trim().split(/\s+/).join('\t');
        }
        add('--num-fragments', cc_numFragments);
        add('--simulate-mode', cc_simulateMode);
        // Distractor mode: CT or fraction
        if (cc_distractorMode === 'ct') {
          if (cc_ctSweep && cc_ctList.trim()) {
            a['--ct-values'] = cc_ctList.trim().split(/\s+/).join('\t');
          } else if (cc_ctFixed) {
            add('--ct', cc_ctFixed);
          }
          // CT calibration
          if (cc_ctCalibMode === 'two-point') {
            a['--ct-calibration'] = `${cc_ctCal1}\t${cc_ctCal2}`;
          } else {
            if (cc_ctBaseline !== '20') add('--ct-baseline', cc_ctBaseline);
            if (cc_ctBaselineFraction !== '0.01') add('--ct-baseline-fraction', cc_ctBaselineFraction);
            if (cc_ctEfficiency !== '1.0') add('--ct-efficiency', cc_ctEfficiency);
          }
        } else {
          if (cc_distractorFracSweep && cc_distractorFracList.trim()) {
            a['--distractor-fraction-values'] = cc_distractorFracList.trim().split(/\s+/).join('\t');
          } else if (cc_distractorFracFixed) {
            add('--distractor-fraction', cc_distractorFracFixed);
          }
        }
        // Hybridization temperature (thermodynamic only)
        if (cc_simulateMode === 'thermodynamic') {
          if (cc_tempSweep && cc_tempList.trim()) {
            a['--hybridization-temperature-values'] = cc_tempList.trim().split(/\s+/).join('\t');
          } else if (cc_tempFixed) {
            add('--hybridization-temperature', cc_tempFixed);
          }
        }
        // Capture fraction
        if (cc_cfSweep && cc_cfList.trim()) {
          a['--capture-fraction-values'] = cc_cfList.trim().split(/\s+/).join('\t');
        } else {
          add('--capture-fraction', cc_cfFixed);
        }
        // Num sequences
        if (cc_nsSweep && cc_nsList.trim()) {
          a['--num-sequences-values'] = cc_nsList.trim().split(/\s+/).join('\t');
        } else if (cc_nsFixed) {
          add('--num-sequences', cc_nsFixed);
        }
        if (cc_seed) add('--seed', cc_seed);
        if (cc_hostFasta) add('--host-fasta', cc_hostFasta);
        if (cc_minimapPreset !== 'sr') add('--minimap-preset', cc_minimapPreset);
        if (cc_hostMinimapPreset !== 'sr') add('--host-minimap-preset', cc_hostMinimapPreset);
        if (cc_fragLenMean !== '175') add('--fragment-length-mean', cc_fragLenMean);
        if (cc_fragLenMin !== '150') add('--fragment-length-min', cc_fragLenMin);
        if (cc_fragLenMax !== '200') add('--fragment-length-max', cc_fragLenMax);
        add('--threads', threads);
        break;
      }
    }

    return a;
  }

  async function runTool() {
    runError = '';
    const args = buildArgs();
    if (!args) return;

    // Persist form state before launching
    try {
      const store = await getFormStore();
      await store.set(`formState_${selectedTool}`, captureState(selectedTool));
      await store.set('lastTool', selectedTool);
      await store.save();
    } catch { /* non-fatal */ }

    const config: PipelineConfig = {
      tool: selectedTool,
      args,
      conda_env: $condaEnvPath,
    };

    launching = true;
    logLines.set([]);
    reportPath.set(null);
    pipelineStatus.set('running');

    // Infer report path from outdir for the "run" tool
    if (outdir && selectedTool === 'run') {
      reportPath.set(`${outdir}/report.html`);
    }

    try {
      await invoke('run_pipeline', { config });
      currentView.set('log');
    } catch (e) {
      runError = String(e);
      pipelineStatus.set('idle');
    } finally {
      launching = false;
    }
  }

  function goToSetup() {
    currentView.set('setup');
  }

  $: isRunning = $pipelineStatus === 'running';
</script>

<div class="run-view">
  <!-- Header -->
  <div class="header">
    <span class="logo">BaitBench</span>
    <div class="spacer"></div>
    <button class="btn-ghost small" on:click={() => open('https://niel-infante.github.io/BaitBench/')}>Documentation</button>
    <button class="btn-ghost small" on:click={goToSetup}>⚙ Change Environment</button>
    {#if isRunning}
      <button class="btn-ghost small active-run" on:click={() => currentView.set('log')}>
        ● View Running Pipeline
      </button>
    {/if}
  </div>

  <div class="body">
    <!-- Tool selector sidebar -->
    <nav class="sidebar">
      {#each CATEGORIES as cat}
        <div class="cat-label">{cat}</div>
        {#each TOOLS.filter(t => t.category === cat) as t}
          <button
            class="tool-btn"
            class:active={selectedTool === t.id}
            on:click={() => switchTool(t.id)}
          >{t.label}</button>
        {/each}
      {/each}
    </nav>

    <!-- Form area -->
    <div class="form-area">
      <div class="tool-header">
        <h2>{tool.label}</h2>
        <p class="tool-desc">{tool.description}</p>
      </div>

      <div class="form-scroll">
        <!-- ── run ─────────────────────────────────────── -->
        {#if selectedTool === 'run'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={r_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={r_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Distractor FASTA(s)" bind:value={r_distractors} multiple required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>

          <section class="form-section">
            <h3>Sample Manifest</h3>
            <div class="toggle-row">
              <button class="toggle-btn" class:active={r_sampleIsFile}
                on:click={() => r_sampleIsFile = true}>File (TSV)</button>
              <button class="toggle-btn" class:active={!r_sampleIsFile}
                on:click={() => r_sampleIsFile = false}>Inline IDs</button>
            </div>
            {#if r_sampleIsFile}
              <FilePicker label="Sample manifest" bind:value={r_sample}
                filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            {:else}
              <label class="field-label" for="sample-inline">
                Space-separated target IDs (optionally followed by weight)
              </label>
              <input id="sample-inline" class="text-input" type="text"
                bind:value={r_sampleInline}
                placeholder="target_1 target_2 3 target_3" />
            {/if}
          </section>

          <section class="form-section">
            <h3>Distractor Fraction</h3>
            <div class="toggle-row">
              <button class="toggle-btn" class:active={r_distractorMode === 'fraction'}
                on:click={() => r_distractorMode = 'fraction'}>Fraction</button>
              <button class="toggle-btn" class:active={r_distractorMode === 'ct'}
                on:click={() => r_distractorMode = 'ct'}>CT Value</button>
            </div>
            {#if r_distractorMode === 'fraction'}
              <div class="field-row">
                <label class="field-label" for="dfrac">Distractor fraction (0–1)</label>
                <input id="dfrac" class="text-input short" type="number" min="0" max="1" step="0.01"
                  bind:value={r_distractorFraction} />
              </div>
            {:else}
              <div class="warning-callout">
                CT scoring is experimental. Results may not reflect real-world assay performance.
              </div>
              <div class="field-group">
                <div class="field-row">
                  <label class="field-label" for="ct">CT value</label>
                  <input id="ct" class="text-input short" type="number" step="0.1"
                    bind:value={r_ct} />
                </div>
                <AdvancedOptions label="CT calibration">
                  <div class="toggle-row" style="margin-bottom:8px">
                    <button class="toggle-btn" class:active={r_ctCalibMode === 'single'}
                      on:click={() => r_ctCalibMode = 'single'}>Single-point</button>
                    <button class="toggle-btn" class:active={r_ctCalibMode === 'two-point'}
                      on:click={() => r_ctCalibMode = 'two-point'}>Two-point</button>
                  </div>
                  {#if r_ctCalibMode === 'single'}
                    <div class="field-row">
                      <label class="field-label" for="ctb">CT baseline</label>
                      <input id="ctb" class="text-input short" type="number" step="0.1"
                        bind:value={r_ctBaseline} />
                    </div>
                    <div class="field-row">
                      <label class="field-label" for="ctbf">CT baseline fraction</label>
                      <input id="ctbf" class="text-input short" type="number" step="0.001" min="0" max="1"
                        bind:value={r_ctBaselineFraction} />
                    </div>
                    <div class="field-row">
                      <label class="field-label" for="cte">PCR efficiency (0–1)</label>
                      <input id="cte" class="text-input short" type="number" step="0.01" min="0" max="1"
                        bind:value={r_ctEfficiency} />
                    </div>
                  {:else}
                    <p class="hint-sm">Provide two (CT, target-fraction) reference points; efficiency is derived automatically.</p>
                    <div class="field-row">
                      <label class="field-label" for="ctcal1">Point 1 (CT,fraction)</label>
                      <input id="ctcal1" class="text-input short" type="text"
                        bind:value={r_ctCal1} placeholder="20,0.01" />
                    </div>
                    <div class="field-row">
                      <label class="field-label" for="ctcal2">Point 2 (CT,fraction)</label>
                      <input id="ctcal2" class="text-input short" type="text"
                        bind:value={r_ctCal2} placeholder="25,0.003" />
                    </div>
                  {/if}
                </AdvancedOptions>
              </div>
            {/if}
          </section>

          <section class="form-section">
            <h3>Simulation</h3>
            <div class="field-row">
              <label class="field-label" for="simmode">Simulate mode</label>
              <select id="simmode" class="select-input" bind:value={r_simulateMode}>
                <option value="thermodynamic">Thermodynamic (TNN)</option>
                <option value="simple">Simple</option>
              </select>
            </div>
            {#if r_simulateMode === 'thermodynamic'}
              <div class="field-row">
                <label class="field-label" for="hybtemp">Hybridization temperature (°C)</label>
                <input id="hybtemp" class="text-input short" type="number" step="1"
                  bind:value={r_hybTemp} />
              </div>
            {/if}
            <div class="field-row">
              <label class="field-label" for="nfrags">Number of fragments</label>
              <input id="nfrags" class="text-input short" type="number" min="100"
                bind:value={r_numFragments} />
            </div>
            <div class="field-row">
              <label class="field-label" for="capfrac">Capture fraction</label>
              <input id="capfrac" class="text-input short" type="number" min="0" max="1" step="0.01"
                bind:value={r_captureFraction} />
            </div>
            <div class="field-row">
              <label class="field-label" for="readlen">Read length (bp)</label>
              <input id="readlen" class="text-input short" type="number" min="1"
                bind:value={r_readLength} />
            </div>
          </section>

          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="report-mode">Report</label>
              <select id="report-mode" class="select-input" bind:value={reportMode}>
                <option value="full">Full HTML</option>
                <option value="both-r">Full HTML + RMarkdown source</option>
                <option value="rmd">RMarkdown source only</option>
                <option value="none">None</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="r-outprefix">Output prefix</label>
              <input id="r-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={cleanup} />
              Clean up intermediate files after run
            </label>
          </section>

          <AdvancedOptions label="Advanced Options">
            <div class="field-row">
              <label class="field-label" for="threads">Threads</label>
              <input id="threads" class="text-input short" type="number" min="1"
                bind:value={threads} />
            </div>
            <div class="field-row">
              <label class="field-label" for="seed">Random seed (blank = random)</label>
              <input id="seed" class="text-input short" type="text" bind:value={r_seed}
                placeholder="e.g. 42" />
            </div>
            <div class="field-row">
              <label class="field-label" for="runname">Run name</label>
              <input id="runname" class="text-input" type="text" bind:value={r_runName} />
            </div>
            <div class="field-row">
              <label class="field-label" for="fraglenmean">Fragment length mean (bp)</label>
              <input id="fraglenmean" class="text-input short" type="number" bind:value={r_fragLenMean} />
            </div>
            <div class="field-row">
              <label class="field-label" for="fraglenmin">Fragment length min (bp)</label>
              <input id="fraglenmin" class="text-input short" type="number" bind:value={r_fragLenMin} />
            </div>
            <div class="field-row">
              <label class="field-label" for="fraglenmax">Fragment length max (bp)</label>
              <input id="fraglenmax" class="text-input short" type="number" bind:value={r_fragLenMax} />
            </div>
            <div class="field-row">
              <label class="field-label" for="minimap">Minimap2 preset</label>
              <input id="minimap" class="text-input short" type="text" bind:value={r_minimapPreset} />
            </div>
            <h4 class="subsection">Sequencing</h4>
            <div class="field-row">
              <label class="field-label" for="r-numseq">Num sequences to sample (blank = all)</label>
              <input id="r-numseq" class="text-input short" type="number" min="1"
                bind:value={r_numSequences} placeholder="all" />
            </div>
            <div class="field-row">
              <label class="field-label" for="r-outfmt">Output format</label>
              <select id="r-outfmt" class="select-input" bind:value={r_outputFormat}>
                <option value="fasta">FASTA</option>
                <option value="fastq">FASTQ</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="r-readsim">Read simulator</label>
              <select id="r-readsim" class="select-input" bind:value={r_readSimulator}>
                <option value="perfect">Perfect (no errors)</option>
                <option value="art">ART (Illumina profiles)</option>
                <option value="badread">Badread (long reads)</option>
              </select>
            </div>
            {#if r_readSimulator === 'art'}
              <div class="field-row">
                <label class="field-label" for="r-artprofile">Sequencer profile</label>
                <input id="r-artprofile" class="text-input" type="text" bind:value={r_sequencerProfile}
                  placeholder="e.g. HS25" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-covdepth">Coverage depth</label>
                <input id="r-covdepth" class="text-input short" type="number" step="0.1" min="0"
                  bind:value={r_coverageDepth} />
              </div>
              <label class="check-label">
                <input type="checkbox" bind:checked={r_pairedEnd} />
                Paired-end output
              </label>
              {#if r_pairedEnd}
                <div class="field-row">
                  <label class="field-label" for="r-pefragmean">Insert size mean (bp)</label>
                  <input id="r-pefragmean" class="text-input short" type="number"
                    bind:value={r_peFargLenMean} />
                </div>
                <div class="field-row">
                  <label class="field-label" for="r-pefragsd">Insert size SD (bp)</label>
                  <input id="r-pefragsd" class="text-input short" type="number"
                    bind:value={r_peFargLenSd} />
                </div>
              {/if}
            {:else if r_readSimulator === 'badread'}
              <div class="field-row">
                <label class="field-label" for="r-brprofile">Sequencer profile / model</label>
                <input id="r-brprofile" class="text-input" type="text" bind:value={r_sequencerProfile}
                  placeholder="e.g. nanopore2023" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brcovdepth">Coverage depth</label>
                <input id="r-brcovdepth" class="text-input short" type="number" step="0.1" min="0"
                  bind:value={r_coverageDepth} />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brlenm">Read length mean (bp)</label>
                <input id="r-brlenm" class="text-input short" type="number" min="1"
                  bind:value={r_longReadLenMean} placeholder="default" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brlens">Read length SD (bp)</label>
                <input id="r-brlens" class="text-input short" type="number" min="0"
                  bind:value={r_longReadLenSd} placeholder="default" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brglitch">Glitches (rate,size,skips)</label>
                <input id="r-brglitch" class="text-input short" type="text"
                  bind:value={r_badreadGlitches} placeholder="e.g. 10000,25,5" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brjunk">Junk reads (%)</label>
                <input id="r-brjunk" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_badreadJunkReads} placeholder="0" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brrandom">Random reads (%)</label>
                <input id="r-brrandom" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_badreadRandomReads} placeholder="0" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brchim">Chimeric reads (%)</label>
                <input id="r-brchim" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_badreadChimeras} placeholder="0" />
              </div>
            {/if}
            <h4 class="subsection">Genome Mode (optional)</h4>
            <FilePicker label="Genomes FASTA" bind:value={r_genomes}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map" bind:value={r_sampleTargetMap}
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <h4 class="subsection">Grouping (optional)</h4>
            <FilePicker label="Target groups" bind:value={r_groups}
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <FilePicker label="Distractor groups" bind:value={r_distractorGroups}
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <h4 class="subsection">Host Filtering (optional)</h4>
            <FilePicker label="Host FASTA" bind:value={r_hostFasta}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <h4 class="subsection">Species Identification</h4>
            <label class="check-label">
              <input type="checkbox" bind:checked={r_identify} />
              Run species-level identification after metrics (genome mode)
            </label>
            {#if r_identify}
              <div class="field-row">
                <label class="field-label" for="r-identthresh">Identity threshold (%)</label>
                <input id="r-identthresh" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_identifyIdentThreshold} />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-minuniq">Min unique targets for call</label>
                <input id="r-minuniq" class="text-input short" type="number" min="1"
                  bind:value={r_identifyMinUniq} />
              </div>
            {/if}
          </AdvancedOptions>

        <!-- ── build-probes ──────────────────────────── -->
        {:else if selectedTool === 'build-probes'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={bp_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section">
            <h3>Design Method</h3>
            <div class="field-row">
              <label class="field-label" for="bp-method">Method</label>
              <select id="bp-method" class="select-input" bind:value={bp_method}>
                <option value="tile">Tile (sliding window)</option>
                <option value="catch-lite">CATCH-lite (native Rust)</option>
                <option value="syotti-lite">Syotti-lite (native Rust)</option>
                <option value="probetools-lite">ProbeTools-lite (native Rust)</option>
                <option value="catch">CATCH (external, requires conda pkg)</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-probelen">Probe length (bp)</label>
              <input id="bp-probelen" class="text-input short" type="number" min="1"
                bind:value={bp_probeLength} />
            </div>
            {#if bp_method === 'tile'}
              <div class="field-row">
                <label class="field-label" for="bp-step">Step (negative = overlap)</label>
                <input id="bp-step" class="text-input short" type="number" bind:value={bp_step} />
              </div>
            {:else if bp_method === 'catch-lite' || bp_method === 'catch'}
              <div class="field-row">
                <label class="field-label" for="bp-catchstride">Probe stride</label>
                <input id="bp-catchstride" class="text-input short" type="number" bind:value={bp_catchStride} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchmm">Mismatches</label>
                <input id="bp-catchmm" class="text-input short" type="number" bind:value={bp_catchMismatches} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchext">Extension (bp each side)</label>
                <input id="bp-catchext" class="text-input short" type="number" min="0" bind:value={bp_catchExtension} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchcov">Min coverage fraction</label>
                <input id="bp-catchcov" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_catchCoverage} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchmh">MinHash dedup threshold</label>
                <input id="bp-catchmh" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_catchMinhashThreshold} />
              </div>
            {:else if bp_method === 'syotti-lite'}
              <div class="field-row">
                <label class="field-label" for="bp-syottimm">Mismatches</label>
                <input id="bp-syottimm" class="text-input short" type="number" bind:value={bp_syottiMismatches} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-syottiseed">Seed length</label>
                <input id="bp-syottiseed" class="text-input short" type="number" bind:value={bp_syottiSeedLen} />
              </div>
            {:else if bp_method === 'probetools-lite'}
              <div class="field-row">
                <label class="field-label" for="bp-ptstep">K-mer step</label>
                <input id="bp-ptstep" class="text-input short" type="number" min="1" bind:value={bp_ptStep} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptident">Cluster identity threshold</label>
                <input id="bp-ptident" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_ptIdentity} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptcov">Target coverage fraction</label>
                <input id="bp-ptcov" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_ptCoverage} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptbatch">Probes per iteration</label>
                <input id="bp-ptbatch" class="text-input short" type="number" min="1" bind:value={bp_ptBatchSize} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptmax">Max panel size (blank = unlimited)</label>
                <input id="bp-ptmax" class="text-input short" type="number" min="1" bind:value={bp_ptMaxPanelSize} placeholder="unlimited" />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptdepth">Min coverage depth</label>
                <input id="bp-ptdepth" class="text-input short" type="number" min="1" bind:value={bp_ptMinDepth} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptiter">Max iterations</label>
                <input id="bp-ptiter" class="text-input short" type="number" min="1" bind:value={bp_ptMaxIterations} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptgain">Min coverage gain to continue</label>
                <input id="bp-ptgain" class="text-input short" type="number" step="0.0001" min="0" bind:value={bp_ptMinCoverageGain} />
              </div>
            {/if}
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="bp-report">Report</label>
              <select id="bp-report" class="select-input" bind:value={reportMode}>
                <option value="full">Full HTML</option>
                <option value="both-r">Full HTML + RMarkdown source</option>
                <option value="rmd">RMarkdown source only</option>
                <option value="none">None</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-outprefix">Output prefix</label>
              <input id="bp-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={cleanup} />
              Clean up intermediate files after run
            </label>
          </section>
          <AdvancedOptions label="Filtering & Assessment">
            <div class="field-row">
              <label class="field-label" for="bp-mingc">Min GC</label>
              <input id="bp-mingc" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_minGc} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-maxgc">Max GC</label>
              <input id="bp-maxgc" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_maxGc} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-maxnfrac">Max N fraction in targets</label>
              <input id="bp-maxnfrac" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_maxNFrac} />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={bp_noNInProbes} />
              Replace N bases in probes with non-N
            </label>
            <div class="field-row">
              <label class="field-label" for="bp-dust">sDUST threshold</label>
              <input id="bp-dust" class="text-input short" type="number" step="0.1" bind:value={bp_dustThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-dustwin">sDUST window size (bp)</label>
              <input id="bp-dustwin" class="text-input short" type="number" min="1" bind:value={bp_dustWindow} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-maxmask">Max masked fraction</label>
              <input id="bp-maxmask" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_maxMaskedFrac} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-collapse">Collapse identity threshold</label>
              <input id="bp-collapse" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_collapseThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-dedup">Dedup identity threshold</label>
              <input id="bp-dedup" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_dedupThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-thresh">Xreact homology threshold (%)</label>
              <input id="bp-thresh" class="text-input short" type="number" min="0" max="100" bind:value={bp_threshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-prox">Proximity distance (bp)</label>
              <input id="bp-prox" class="text-input short" type="number" min="0" bind:value={bp_proximity} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-minimap">Minimap2 preset</label>
              <input id="bp-minimap" class="text-input short" type="text" bind:value={bp_minimapPreset} />
            </div>
            <h4 class="subsection">Coverage Refinement</h4>
            <div class="field-row">
              <label class="field-label" for="bp-refmode">Refinement</label>
              <select id="bp-refmode" class="select-input" bind:value={bp_refineMode}>
                <option value="none">None</option>
                <option value="iterations">Fixed iterations</option>
                <option value="stable">Until stable</option>
              </select>
            </div>
            {#if bp_refineMode !== 'none'}
              <div class="field-row">
                <label class="field-label" for="bp-refthresh">1× coverage threshold (%)</label>
                <input id="bp-refthresh" class="text-input short" type="number" min="0" max="100" bind:value={bp_refineThreshold} />
              </div>
              {#if bp_refineMode === 'iterations'}
                <div class="field-row">
                  <label class="field-label" for="bp-refiter">Iterations</label>
                  <input id="bp-refiter" class="text-input short" type="number" min="1" bind:value={bp_refineIterations} />
                </div>
              {/if}
            {/if}
            <h4 class="subsection">Assessment</h4>
            <FilePicker label="Genomes FASTA (cross-reactivity check)" bind:value={bp_genomes}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <label class="check-label">
              <input type="checkbox" bind:checked={bp_skipAssess} />
              Skip assessment step
            </label>
            <div class="field-row">
              <label class="field-label" for="bp-threads">Threads</label>
              <input id="bp-threads" class="text-input short" type="number" min="1" bind:value={threads} />
            </div>
          </AdvancedOptions>

        <!-- ── assess-probes ────────────────────────── -->
        {:else if selectedTool === 'assess-probes'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={ap_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={ap_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Genomes FASTA (optional)" bind:value={ap_genomes}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="ap-report">Report</label>
              <select id="ap-report" class="select-input" bind:value={reportMode}>
                <option value="full">Full HTML</option>
                <option value="both-r">Full HTML + RMarkdown source</option>
                <option value="rmd">RMarkdown source only</option>
                <option value="none">None</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="ap-outprefix">Output prefix</label>
              <input id="ap-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={cleanup} />
              Clean up intermediate files after run
            </label>
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="ap-threshold">Homology threshold (%)</label>
              <input id="ap-threshold" class="text-input short" type="number" min="0" max="100"
                bind:value={ap_threshold} />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={ap_noIndividual} />
              Skip per-target individual coverage mapping (faster for large panels)
            </label>
            <div class="field-row">
              <label class="field-label" for="ap-prox">Proximity distance (bp)</label>
              <input id="ap-prox" class="text-input short" type="number" min="0" bind:value={ap_proximity} />
            </div>
            <div class="field-row">
              <label class="field-label" for="ap-minimap">Minimap2 preset</label>
              <input id="ap-minimap" class="text-input short" type="text" bind:value={ap_minimapPreset} />
            </div>
            <div class="field-row">
              <label class="field-label" for="ap-gaplen">Min gap length in detail output (bp)</label>
              <input id="ap-gaplen" class="text-input short" type="number" min="1"
                bind:value={ap_gapMinLength} placeholder="default" />
            </div>
            <h4 class="subsection">Coverage Refinement</h4>
            <div class="field-row">
              <label class="field-label" for="ap-refmode">Refinement</label>
              <select id="ap-refmode" class="select-input" bind:value={ap_refineMode}>
                <option value="none">None</option>
                <option value="iterations">Fixed iterations</option>
                <option value="stable">Until stable</option>
              </select>
            </div>
            {#if ap_refineMode !== 'none'}
              <div class="field-row">
                <label class="field-label" for="ap-refthresh">1× coverage threshold (%)</label>
                <input id="ap-refthresh" class="text-input short" type="number" min="0" max="100" bind:value={ap_refineThreshold} />
              </div>
              {#if ap_refineMode === 'iterations'}
                <div class="field-row">
                  <label class="field-label" for="ap-refiter">Iterations</label>
                  <input id="ap-refiter" class="text-input short" type="number" min="1" bind:value={ap_refineIterations} />
                </div>
              {/if}
            {/if}
            <div class="field-row">
              <label class="field-label" for="ap-threads">Threads</label>
              <input id="ap-threads" class="text-input short" type="number" min="1" bind:value={threads} />
            </div>
          </AdvancedOptions>

        <!-- ── xreact ───────────────────────────────── -->
        {:else if selectedTool === 'xreact'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Probes FASTA" bind:value={xr_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Against FASTA(s) (genomes to check against)" multiple
              bind:value={xr_against}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <label class="check-label">
              <input type="checkbox" bind:checked={xr_self} />
              Probe-vs-probe self cross-reactivity
            </label>
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="xr-report">Report</label>
              <select id="xr-report" class="select-input" bind:value={reportMode}>
                <option value="full">Full HTML</option>
                <option value="both-r">Full HTML + RMarkdown source</option>
                <option value="rmd">RMarkdown source only</option>
                <option value="none">None</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="xr-outprefix">Output prefix</label>
              <input id="xr-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={cleanup} />
              Clean up intermediate files after run
            </label>
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="xr-threshold">Homology threshold (%)</label>
              <input id="xr-threshold" class="text-input short" type="number" min="0" max="100"
                bind:value={xr_threshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="xr-minimap">Minimap2 preset</label>
              <input id="xr-minimap" class="text-input short" type="text" bind:value={xr_minimapPreset} />
            </div>
          </AdvancedOptions>

        <!-- ── panel-qc ─────────────────────────────── -->
        {:else if selectedTool === 'panel-qc'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={pq_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map (TSV)" bind:value={pq_sampleTargetMap} required
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="pq-report">Report</label>
              <select id="pq-report" class="select-input" bind:value={reportMode}>
                <option value="full">Full HTML</option>
                <option value="both-r">Full HTML + RMarkdown source</option>
                <option value="rmd">RMarkdown source only</option>
                <option value="none">None</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="pq-outprefix">Output prefix</label>
              <input id="pq-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={cleanup} />
              Clean up intermediate files after run
            </label>
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="pq-ident">Identity threshold (%)</label>
              <input id="pq-ident" class="text-input short" type="number" min="0" max="100"
                bind:value={pq_identityThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="pq-minimap">Minimap2 preset</label>
              <input id="pq-minimap" class="text-input short" type="text" bind:value={pq_minimapPreset} />
            </div>
          </AdvancedOptions>

        <!-- ── identify ─────────────────────────────── -->
        {:else if selectedTool === 'identify'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="detected_detail.tsv" bind:value={id_detectedDetail} required
              filters={[{ name: 'TSV', extensions: ['tsv'] }]} />
            <FilePicker label="Sample-target map (TSV)" bind:value={id_sampleTargetMap} required
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <FilePicker label="Pre-computed target similarity TSV (optional)"
              bind:value={id_targetSimilarity}
              filters={[{ name: 'TSV', extensions: ['tsv'] }]} />
            <FilePicker label="Targets FASTA (optional, compute similarity on-the-fly)"
              bind:value={id_targets}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="id-outprefix">Output prefix</label>
              <input id="id-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="id-ident">Identity threshold (%)</label>
              <input id="id-ident" class="text-input short" type="number" min="0" max="100"
                bind:value={id_identityThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="id-minuniq">Min unique targets for call</label>
              <input id="id-minuniq" class="text-input short" type="number" min="1"
                bind:value={id_minUniqueTargets} />
            </div>
            <div class="field-row">
              <label class="field-label" for="id-minimap">Minimap2 preset</label>
              <input id="id-minimap" class="text-input short" type="text" bind:value={id_minimapPreset} />
            </div>
          </AdvancedOptions>

        <!-- ── coverage-curve ────────────────────────── -->
        {:else if selectedTool === 'coverage-curve'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={cc_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={cc_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Distractor FASTA(s)" bind:value={cc_distractors} multiple required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Genomes FASTA (optional, genome mode)" bind:value={cc_genomes}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map (optional)" bind:value={cc_sampleTargetMap}
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
          </section>

          <section class="form-section">
            <h3>Sample Manifest</h3>
            <div class="toggle-row">
              <button class="toggle-btn" class:active={cc_sampleIsFile}
                on:click={() => cc_sampleIsFile = true}>File (TSV)</button>
              <button class="toggle-btn" class:active={!cc_sampleIsFile}
                on:click={() => cc_sampleIsFile = false}>Inline IDs</button>
            </div>
            {#if cc_sampleIsFile}
              <FilePicker label="Sample manifest (optional)" bind:value={cc_sample}
                filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            {:else}
              <label class="field-label" for="cc-sample-inline">
                Space-separated target IDs (optionally followed by weight)
              </label>
              <input id="cc-sample-inline" class="text-input" type="text"
                bind:value={cc_sampleInline}
                placeholder="target_1 target_2 target_3" />
            {/if}
          </section>

          <section class="form-section">
            <h3>Simulation</h3>
            <div class="field-row">
              <label class="field-label" for="cc-nfrags">Number of fragments</label>
              <input id="cc-nfrags" class="text-input short" type="number" min="100"
                bind:value={cc_numFragments} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-simmode">Simulate mode</label>
              <select id="cc-simmode" class="select-input" bind:value={cc_simulateMode}>
                <option value="thermodynamic">Thermodynamic (TNN)</option>
                <option value="simple">Simple</option>
              </select>
            </div>
          </section>

          <section class="form-section">
            <h3>Parameter Sweep</h3>
            <p class="sweep-hint">Check a parameter to sweep over multiple values; uncheck to use a single fixed value.</p>

            <div class="field-row" style="margin-bottom:8px">
              <label class="field-label">Distractor source</label>
              <div class="toggle-row">
                <button class="toggle-btn" class:active={cc_distractorMode === 'ct'}
                  on:click={() => cc_distractorMode = 'ct'}>CT value</button>
                <button class="toggle-btn" class:active={cc_distractorMode === 'fraction'}
                  on:click={() => cc_distractorMode = 'fraction'}>Distractor fraction</button>
              </div>
            </div>

            {#if cc_distractorMode === 'ct'}
              <div class="warning-callout" style="margin-bottom:8px">
                CT scoring is experimental. Results may not reflect real-world assay performance.
              </div>
              <div class="sweep-row">
                <label class="sweep-label">
                  <input type="checkbox" bind:checked={cc_ctSweep} />
                  CT
                </label>
                {#if cc_ctSweep}
                  <input class="text-input sweep-input" type="text" bind:value={cc_ctList} placeholder="20 25 30" />
                {:else}
                  <input class="text-input short" type="number" step="0.1" bind:value={cc_ctFixed} placeholder="25" />
                {/if}
              </div>
              <AdvancedOptions label="CT calibration">
                <div class="toggle-row" style="margin-bottom:8px">
                  <button class="toggle-btn" class:active={cc_ctCalibMode === 'single'}
                    on:click={() => cc_ctCalibMode = 'single'}>Single-point</button>
                  <button class="toggle-btn" class:active={cc_ctCalibMode === 'two-point'}
                    on:click={() => cc_ctCalibMode = 'two-point'}>Two-point</button>
                </div>
                {#if cc_ctCalibMode === 'single'}
                  <div class="field-row">
                    <label class="field-label" for="cc-ctb">CT baseline</label>
                    <input id="cc-ctb" class="text-input short" type="number" step="0.1" bind:value={cc_ctBaseline} />
                  </div>
                  <div class="field-row">
                    <label class="field-label" for="cc-ctbf">CT baseline fraction</label>
                    <input id="cc-ctbf" class="text-input short" type="number" step="0.001" min="0" max="1" bind:value={cc_ctBaselineFraction} />
                  </div>
                  <div class="field-row">
                    <label class="field-label" for="cc-cte">PCR efficiency (0–1)</label>
                    <input id="cc-cte" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={cc_ctEfficiency} />
                  </div>
                {:else}
                  <p class="hint-sm">Provide two (CT, target-fraction) reference points; efficiency is derived automatically.</p>
                  <div class="field-row">
                    <label class="field-label" for="cc-ctcal1">Point 1 (CT,fraction)</label>
                    <input id="cc-ctcal1" class="text-input short" type="text" bind:value={cc_ctCal1} placeholder="20,0.01" />
                  </div>
                  <div class="field-row">
                    <label class="field-label" for="cc-ctcal2">Point 2 (CT,fraction)</label>
                    <input id="cc-ctcal2" class="text-input short" type="text" bind:value={cc_ctCal2} placeholder="25,0.003" />
                  </div>
                {/if}
              </AdvancedOptions>
            {:else}
              <div class="sweep-row">
                <label class="sweep-label">
                  <input type="checkbox" bind:checked={cc_distractorFracSweep} />
                  Distractor Fraction
                </label>
                {#if cc_distractorFracSweep}
                  <input class="text-input sweep-input" type="text" bind:value={cc_distractorFracList} placeholder="0.5 0.7 0.9" />
                {:else}
                  <input class="text-input short" type="number" min="0" max="1" step="0.01" bind:value={cc_distractorFracFixed} placeholder="0.9" />
                {/if}
              </div>
            {/if}

            {#if cc_simulateMode === 'thermodynamic'}
              <div class="sweep-row">
                <label class="sweep-label">
                  <input type="checkbox" bind:checked={cc_tempSweep} />
                  Hybridization Temp (°C)
                </label>
                {#if cc_tempSweep}
                  <input class="text-input sweep-input" type="text" bind:value={cc_tempList} placeholder="60 65 70 75" />
                {:else}
                  <input class="text-input short" type="number" step="1" bind:value={cc_tempFixed} placeholder="70" />
                {/if}
              </div>
            {/if}

            <div class="sweep-row">
              <label class="sweep-label">
                <input type="checkbox" bind:checked={cc_cfSweep} />
                Capture Fraction
              </label>
              {#if cc_cfSweep}
                <input class="text-input sweep-input" type="text" bind:value={cc_cfList} placeholder="0.3 0.5 0.8" />
              {:else}
                <input class="text-input short" type="number" min="0" max="1" step="0.01" bind:value={cc_cfFixed} placeholder="0.5" />
              {/if}
            </div>

            <div class="sweep-row">
              <label class="sweep-label">
                <input type="checkbox" bind:checked={cc_nsSweep} />
                Num Sequences
              </label>
              {#if cc_nsSweep}
                <input class="text-input sweep-input" type="text" bind:value={cc_nsList} placeholder="500 1000 5000" />
              {:else}
                <input class="text-input short" type="number" min="1" bind:value={cc_nsFixed} placeholder="all" />
              {/if}
            </div>
          </section>

          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="cc-report">Report</label>
              <select id="cc-report" class="select-input" bind:value={reportMode}>
                <option value="full">Full HTML</option>
                <option value="both-r">Full HTML + RMarkdown source</option>
                <option value="rmd">RMarkdown source only</option>
                <option value="none">None</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-outprefix">Output prefix</label>
              <input id="cc-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={cleanup} />
              Clean up intermediate files after run
            </label>
          </section>
          <AdvancedOptions label="Advanced">
            <div class="field-row">
              <label class="field-label" for="cc-threads">Threads</label>
              <input id="cc-threads" class="text-input short" type="number" min="1" bind:value={threads} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-seed">Random seed (blank = random)</label>
              <input id="cc-seed" class="text-input short" type="text" bind:value={cc_seed} placeholder="e.g. 42" />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-fraglenmean">Fragment length mean (bp)</label>
              <input id="cc-fraglenmean" class="text-input short" type="number" bind:value={cc_fragLenMean} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-fraglenmin">Fragment length min (bp)</label>
              <input id="cc-fraglenmin" class="text-input short" type="number" bind:value={cc_fragLenMin} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-fraglenmax">Fragment length max (bp)</label>
              <input id="cc-fraglenmax" class="text-input short" type="number" bind:value={cc_fragLenMax} />
            </div>
            <FilePicker label="Host FASTA (optional filtering)" bind:value={cc_hostFasta}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <div class="field-row">
              <label class="field-label" for="cc-minimap">Minimap2 preset</label>
              <input id="cc-minimap" class="text-input short" type="text" bind:value={cc_minimapPreset} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-hostminimap">Host minimap2 preset</label>
              <input id="cc-hostminimap" class="text-input short" type="text" bind:value={cc_hostMinimapPreset} />
            </div>
          </AdvancedOptions>
        {/if}

        <!-- Error -->
        {#if runError}
          <div class="run-error">{runError}</div>
        {/if}

        <!-- Run button -->
        <div class="run-row">
          <button
            class="btn-run"
            disabled={isRunning || launching}
            on:click={runTool}
          >
            {#if launching}
              <span class="spinner-sm"></span> Starting…
            {:else if isRunning}
              Running…
            {:else}
              ▶ Run {tool.label}
            {/if}
          </button>
          {#if isRunning}
            <button class="btn-ghost" on:click={() => currentView.set('log')}>
              View log →
            </button>
          {/if}
        </div>
      </div>
    </div>
  </div>
</div>

<style>
  .run-view {
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
  }
  .header {
    display: flex;
    align-items: center;
    gap: 10px;
    padding: 8px 16px;
    border-bottom: 1px solid var(--color-border);
    background: var(--color-card);
    flex-shrink: 0;
  }
  .logo {
    font-size: 1.05rem;
    font-weight: 800;
    color: var(--color-primary);
  }
  .spacer { flex: 1; }
  .body {
    display: flex;
    flex: 1;
    overflow: hidden;
  }
  /* Sidebar */
  .sidebar {
    width: 180px;
    flex-shrink: 0;
    border-right: 1px solid var(--color-border);
    background: var(--color-card);
    padding: 10px 0;
    overflow-y: auto;
    display: flex;
    flex-direction: column;
    gap: 1px;
  }
  .cat-label {
    font-size: 0.7rem;
    font-weight: 700;
    text-transform: uppercase;
    color: var(--color-muted);
    letter-spacing: 0.04em;
    padding: 8px 12px 3px;
  }
  .tool-btn {
    width: 100%;
    text-align: left;
    padding: 7px 14px;
    font-size: 0.85rem;
    border: none;
    background: transparent;
    cursor: pointer;
    color: var(--color-text);
    border-radius: 0;
  }
  .tool-btn:hover { background: var(--color-btn-hover); }
  .tool-btn.active {
    background: var(--color-primary-light);
    color: var(--color-primary);
    font-weight: 600;
  }
  /* Form area */
  .form-area {
    flex: 1;
    display: flex;
    flex-direction: column;
    overflow: hidden;
  }
  .tool-header {
    padding: 14px 20px 10px;
    border-bottom: 1px solid var(--color-border);
    flex-shrink: 0;
  }
  .tool-header h2 { margin: 0 0 3px; font-size: 1.05rem; }
  .tool-desc { margin: 0; font-size: 0.82rem; color: var(--color-muted); }
  .form-scroll {
    flex: 1;
    overflow-y: auto;
    padding: 16px 20px 32px;
    display: flex;
    flex-direction: column;
    gap: 16px;
  }
  .form-section {
    display: flex;
    flex-direction: column;
    gap: 10px;
  }
  .form-section h3 {
    margin: 0 0 2px;
    font-size: 0.82rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--color-muted);
  }
  .field-row {
    display: flex;
    align-items: center;
    gap: 10px;
  }
  .field-row .field-label {
    flex: 1;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--color-label);
    min-width: 160px;
  }
  .field-group { display: flex; flex-direction: column; gap: 8px; }
  .field-label { font-size: 0.85rem; font-weight: 600; color: var(--color-label); }
  .hint-sm { display: block; font-size: 0.75rem; font-weight: 400; color: var(--color-muted); margin-top: 1px; }
  .sweep-hint { font-size: 0.8rem; color: var(--color-muted); margin: 0 0 8px; }
  .sweep-row { display: flex; align-items: center; gap: 10px; min-height: 32px; }
  .sweep-label { display: flex; align-items: center; gap: 6px; font-size: 0.85rem; font-weight: 600; color: var(--color-label); min-width: 180px; cursor: pointer; white-space: nowrap; }
  .sweep-label input[type="checkbox"] { width: 15px; height: 15px; flex-shrink: 0; cursor: pointer; }
  .sweep-input { flex: 1; min-width: 0; }
  .text-input, .select-input {
    padding: 5px 8px;
    font-size: 0.85rem;
    border: 1px solid var(--color-border);
    border-radius: 4px;
    background: var(--color-input-bg);
    color: var(--color-text);
  }
  .text-input { flex: 1; }
  .text-input.short { width: 90px; flex: none; }
  .select-input { flex: 1; }
  .check-label {
    display: flex;
    align-items: center;
    gap: 8px;
    font-size: 0.85rem;
    color: var(--color-text);
    cursor: pointer;
  }
  .toggle-row {
    display: flex;
    border: 1px solid var(--color-border);
    border-radius: 5px;
    overflow: hidden;
    width: fit-content;
  }
  .toggle-btn {
    padding: 4px 12px;
    font-size: 0.82rem;
    border: none;
    background: var(--color-input-bg);
    cursor: pointer;
    color: var(--color-muted);
  }
  .toggle-btn.active {
    background: var(--color-primary);
    color: #fff;
    font-weight: 600;
  }
  .subsection {
    margin: 8px 0 0;
    font-size: 0.78rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--color-muted);
  }
  .run-error {
    background: #fff5f5;
    border: 1px solid #feb2b2;
    border-radius: 5px;
    padding: 8px 12px;
    font-size: 0.83rem;
    color: #c53030;
  }
  .warning-callout {
    background: #fffbeb;
    border: 1px solid #f6ad55;
    border-radius: 5px;
    padding: 8px 12px;
    font-size: 0.82rem;
    color: #7b341e;
  }
  .run-row {
    display: flex;
    align-items: center;
    gap: 12px;
    margin-top: 4px;
  }
  .btn-run {
    padding: 10px 22px;
    background: var(--color-primary);
    color: #fff;
    border: none;
    border-radius: 6px;
    font-size: 0.92rem;
    font-weight: 700;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 7px;
  }
  .btn-run:disabled { opacity: 0.5; cursor: not-allowed; }
  .btn-run:not(:disabled):hover { opacity: 0.88; }
  .spinner-sm {
    width: 12px; height: 12px;
    border: 2px solid rgba(255,255,255,0.4);
    border-top-color: #fff;
    border-radius: 50%;
    animation: spin 0.7s linear infinite;
  }
  @keyframes spin { to { transform: rotate(360deg); } }
  .btn-ghost {
    padding: 5px 10px;
    background: transparent;
    border: 1px solid var(--color-border);
    border-radius: 5px;
    font-size: 0.82rem;
    cursor: pointer;
    color: var(--color-text);
  }
  .btn-ghost:hover { background: var(--color-btn-hover); }
  .btn-ghost.small { font-size: 0.8rem; padding: 4px 8px; }
  .active-run { color: #2b6cb0; border-color: #bee3f8; }
</style>
