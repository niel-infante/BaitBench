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
  let r_hostMinimapPreset = 'sr';
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
  let bp_aligner = 'minimap2';
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
  let ap_aligner = 'minimap2';
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
  let xr_aligner = 'minimap2';
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
        r_genomes, r_sampleTargetMap, r_groups, r_distractorGroups, r_hostFasta, r_hostMinimapPreset,
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
        bp_proximity, bp_aligner, bp_minimapPreset, bp_threshold,
        bp_refineThreshold, bp_refineMode, bp_refineIterations,
        bp_skipAssess, bp_genomes,
      };
      case 'assess-probes': return { ...shared,
        ap_targets, ap_probes, ap_genomes, ap_threshold, ap_noIndividual,
        ap_proximity, ap_aligner, ap_minimapPreset, ap_gapMinLength,
        ap_refineThreshold, ap_refineMode, ap_refineIterations,
      };
      case 'xreact': return { ...shared,
        xr_probes, xr_against, xr_self, xr_threshold, xr_aligner, xr_minimapPreset,
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
        r_hostMinimapPreset = str('r_hostMinimapPreset', 'sr');
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
        bp_aligner = str('bp_aligner', 'minimap2');
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
        ap_aligner = str('ap_aligner', 'minimap2');
        ap_gapMinLength = str('ap_gapMinLength');
        ap_refineThreshold = str('ap_refineThreshold', '80');
        ap_refineMode = str('ap_refineMode', 'none') as 'none' | 'iterations' | 'stable';
        ap_refineIterations = str('ap_refineIterations', '3');
        break;
      case 'xreact':
        xr_probes = str('xr_probes'); xr_against = str('xr_against');
        xr_self = bl('xr_self'); xr_threshold = str('xr_threshold', '80');
        xr_aligner = str('xr_aligner', 'minimap2');
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
        if (r_hostFasta && r_hostMinimapPreset !== 'sr') add('--host-minimap-preset', r_hostMinimapPreset);
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
        if (bp_aligner !== 'minimap2') add('--aligner', bp_aligner);
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
        if (ap_aligner !== 'minimap2') add('--aligner', ap_aligner);
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
        if (xr_aligner !== 'minimap2') add('--aligner', xr_aligner);
        if (xr_aligner === 'blast') {
          add('--threads', threads);
        } else if (xr_minimapPreset !== 'sr') {
          add('--minimap-preset', xr_minimapPreset);
        }
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
    <div class="header-side"></div>
    <span class="logo">BaitBench</span>
    <div class="header-side right">
      {#if isRunning}
        <button class="btn-running" on:click={() => currentView.set('log')}>
          ● View Running Pipeline
        </button>
      {/if}
    </div>
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
      <div class="sidebar-spacer"></div>
      <div class="sidebar-utilities">
        <button class="util-btn" on:click={() => open('https://niel-infante.github.io/BaitBench/')}>Documentation ↗</button>
        <button class="util-btn" on:click={goToSetup}>⚙ Change Environment</button>
      </div>
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
          <section class="form-section" data-kind="input">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={r_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={r_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Distractor FASTA(s)" bind:value={r_distractors} multiple required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Genomes FASTA (optional, genome mode)" bind:value={r_genomes}
              tooltip="Full genome sequences for genome mode (e.g. complete bacterial chromosomes). Fragments are generated from these genomes but reads are mapped back to Targets FASTA. Use for bacteria or large pathogens where the probe target region differs from the full genome."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map (optional, genome mode)" bind:value={r_sampleTargetMap}
              tooltip="Two-column TSV: genome_id → target_id. Links genome sequences to their probe target regions in genome mode. If omitted, auto-matched by exact name or 'genome_id|target_id' prefix."
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
          </section>

          <section class="form-section" data-kind="input">
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

          <section class="form-section" data-kind="params">
            <h3>Distractor Fraction</h3>
            <div class="toggle-row">
              <button class="toggle-btn" class:active={r_distractorMode === 'fraction'}
                on:click={() => r_distractorMode = 'fraction'}>Fraction</button>
              <button class="toggle-btn" class:active={r_distractorMode === 'ct'}
                on:click={() => r_distractorMode = 'ct'}>CT Value</button>
            </div>
            {#if r_distractorMode === 'fraction'}
              <div class="field-row">
                <label class="field-label" for="dfrac" data-tooltip="Fraction of simulated reads from distractor (background) sequences. 0.9 = 90% background, 10% target. Models real sample composition.">Distractor fraction (0–1) <span class="tip">?</span></label>
                <input id="dfrac" class="text-input short" type="number" min="0" max="1" step="0.01"
                  bind:value={r_distractorFraction} />
              </div>
            {:else}
              <div class="warning-callout">
                CT scoring is experimental. Results may not reflect real-world assay performance.
              </div>
              <div class="field-group">
                <div class="field-row">
                  <label class="field-label" for="ct" data-tooltip="qPCR cycle threshold. Higher CT = more dilute target. CT 20 ≈ 1% target, CT 25 ≈ 0.03%, CT 30 ≈ 0.001% (at 100% PCR efficiency).">CT value <span class="tip">?</span></label>
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
                      <label class="field-label" for="ctb" data-tooltip="CT value of the calibration reference point. Default: CT 20 corresponds to 1% target fraction.">CT baseline <span class="tip">?</span></label>
                      <input id="ctb" class="text-input short" type="number" step="0.1"
                        bind:value={r_ctBaseline} />
                    </div>
                    <div class="field-row">
                      <label class="field-label" for="ctbf" data-tooltip="Target fraction at the baseline CT. Default: 0.01 means CT 20 = 1% target in the mixture.">CT baseline fraction <span class="tip">?</span></label>
                      <input id="ctbf" class="text-input short" type="number" step="0.001" min="0" max="1"
                        bind:value={r_ctBaselineFraction} />
                    </div>
                    <div class="field-row">
                      <label class="field-label" for="cte" data-tooltip="PCR amplification efficiency per cycle. 1.0 = 100% (doubles every cycle). Real assays typically run at 0.90–0.98.">PCR efficiency (0–1) <span class="tip">?</span></label>
                      <input id="cte" class="text-input short" type="number" step="0.01" min="0" max="1"
                        bind:value={r_ctEfficiency} />
                    </div>
                  {:else}
                    <p class="hint-sm">Provide two (CT, target-fraction) reference points; efficiency is derived automatically.</p>
                    <div class="field-row">
                      <label class="field-label" for="ctcal1" data-tooltip="Format: CT,fraction — e.g. '20,0.01' means CT 20 corresponds to 1% target. PCR efficiency is derived from the two points.">Point 1 (CT,fraction) <span class="tip">?</span></label>
                      <input id="ctcal1" class="text-input short" type="text"
                        bind:value={r_ctCal1} placeholder="20,0.01" />
                    </div>
                    <div class="field-row">
                      <label class="field-label" for="ctcal2" data-tooltip="Second calibration point. Format: CT,fraction — e.g. '25,0.003' means CT 25 corresponds to 0.3% target.">Point 2 (CT,fraction) <span class="tip">?</span></label>
                      <input id="ctcal2" class="text-input short" type="text"
                        bind:value={r_ctCal2} placeholder="25,0.003" />
                    </div>
                  {/if}
                </AdvancedOptions>
              </div>
            {/if}
          </section>

          <section class="form-section" data-kind="params">
            <h3>Simulation</h3>
            <div class="field-row">
              <label class="field-label" for="simmode" data-tooltip="Thermodynamic uses the SantaLucia (1998) nearest-neighbor TNN model to score probe-target binding based on hybridization free energy. Simple assigns uniform capture probability. Thermodynamic is more accurate but slower.">Simulate mode <span class="tip">?</span></label>
              <select id="simmode" class="select-input" bind:value={r_simulateMode}>
                <option value="thermodynamic">Thermodynamic (TNN)</option>
                <option value="simple">Simple</option>
              </select>
            </div>
            {#if r_simulateMode === 'thermodynamic'}
              <div class="field-row">
                <label class="field-label" for="hybtemp" data-tooltip="Temperature used in the TNN model to calculate probe-target binding free energy. Higher temp = stricter hybridization. Typical range: 60–75°C.">Hybridization temperature (°C) <span class="tip">?</span></label>
                <input id="hybtemp" class="text-input short" type="number" step="1"
                  bind:value={r_hybTemp} />
              </div>
            {/if}
            <div class="field-row">
              <label class="field-label" for="nfrags" data-tooltip="Total DNA fragments to simulate across all sequences. More fragments → better proportional accuracy but slower run. Typical: 10,000–100,000.">Number of fragments <span class="tip">?</span></label>
              <input id="nfrags" class="text-input short" type="number" min="100"
                bind:value={r_numFragments} />
            </div>
            <div class="field-row">
              <label class="field-label" for="capfrac" data-tooltip="Proportion of fragments overlapping a probe that are retained (captured). 1.0 = capture everything that touches a probe. Models real hybridization capture efficiency.">Capture fraction <span class="tip">?</span></label>
              <input id="capfrac" class="text-input short" type="number" min="0" max="1" step="0.01"
                bind:value={r_captureFraction} />
            </div>
            <div class="field-row">
              <label class="field-label" for="readlen" data-tooltip="Simulated fragments are trimmed to this length to produce reads. Should match the sequencing technology (e.g. 150 bp for Illumina short reads).">Read length (bp) <span class="tip">?</span></label>
              <input id="readlen" class="text-input short" type="number" min="1"
                bind:value={r_readLength} />
            </div>
          </section>

          <section class="form-section" data-kind="output">
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
              <label class="field-label" for="seed" data-tooltip="Fix this value to reproduce the exact same simulation in future runs. Leave blank for a new random seed each time.">Random seed (blank = random) <span class="tip">?</span></label>
              <input id="seed" class="text-input short" type="text" bind:value={r_seed}
                placeholder="e.g. 42" />
            </div>
            <div class="field-row">
              <label class="field-label" for="runname">Run name</label>
              <input id="runname" class="text-input" type="text" bind:value={r_runName} />
            </div>
            <div class="field-row">
              <label class="field-label" for="fraglenmean" data-tooltip="Mean length of simulated DNA fragments before trimming to read length. Longer fragments capture more context per molecule. Default: 500 bp.">Fragment length mean (bp) <span class="tip">?</span></label>
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
              <label class="field-label" for="minimap" data-tooltip="Minimap2 alignment preset for mapping reads back to the reference. 'sr' for Illumina short reads (default), 'map-ont' for Nanopore, 'map-hifi' for PacBio HiFi.">Minimap2 preset <span class="tip">?</span></label>
              <select id="minimap" class="text-input short" bind:value={r_minimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
            </div>
            <h4 class="subsection">Sequencing</h4>
            <div class="field-row">
              <label class="field-label" for="r-numseq" data-tooltip="Subsample this many reads from the simulated pool before mapping. Leave blank to use all reads. Useful for testing at lower sequencing depths.">Num sequences to sample (blank = all) <span class="tip">?</span></label>
              <input id="r-numseq" class="text-input short" type="number" min="1"
                bind:value={r_numSequences} placeholder="all" />
            </div>
            <div class="field-row">
              <label class="field-label" for="r-outfmt" data-tooltip="FASTA for the perfect simulator. FASTQ includes quality scores — required for ART or Badread read simulators.">Output format <span class="tip">?</span></label>
              <select id="r-outfmt" class="select-input" bind:value={r_outputFormat}>
                <option value="fasta">FASTA</option>
                <option value="fastq">FASTQ</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="r-readsim" data-tooltip="Perfect: error-free reads (fastest, default). ART: Illumina error profiles (requires art_modern conda package). Badread: long-read error models for Nanopore/PacBio (requires badread conda package).">Read simulator <span class="tip">?</span></label>
              <select id="r-readsim" class="select-input" bind:value={r_readSimulator}>
                <option value="perfect">Perfect (no errors)</option>
                <option value="art">ART (Illumina profiles)</option>
                <option value="badread">Badread (long reads)</option>
              </select>
            </div>
            {#if r_readSimulator === 'art'}
              <div class="field-row">
                <label class="field-label" for="r-artprofile" data-tooltip="ART instrument profile name. Examples: HS25 (HiSeq 2500), HS10 (HiSeq 1000), NS50 (NextSeq 500). See ART documentation for full list.">Sequencer profile <span class="tip">?</span></label>
                <input id="r-artprofile" class="text-input" type="text" bind:value={r_sequencerProfile}
                  placeholder="e.g. HS25" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-covdepth" data-tooltip="Target sequencing coverage depth (×). ART will generate reads until this average depth is reached per reference sequence.">Coverage depth <span class="tip">?</span></label>
                <input id="r-covdepth" class="text-input short" type="number" step="0.1" min="0"
                  bind:value={r_coverageDepth} />
              </div>
              <label class="check-label">
                <input type="checkbox" bind:checked={r_pairedEnd} />
                Paired-end output
              </label>
              {#if r_pairedEnd}
                <div class="field-row">
                  <label class="field-label" for="r-pefragmean" data-tooltip="Mean DNA fragment insert size for paired-end sequencing (bp). Must be larger than the read length. Typical Illumina: 300–500 bp.">Insert size mean (bp) <span class="tip">?</span></label>
                  <input id="r-pefragmean" class="text-input short" type="number"
                    bind:value={r_peFargLenMean} />
                </div>
                <div class="field-row">
                  <label class="field-label" for="r-pefragsd" data-tooltip="Standard deviation of the fragment insert size (bp). Typical Illumina: 50–100 bp.">Insert size SD (bp) <span class="tip">?</span></label>
                  <input id="r-pefragsd" class="text-input short" type="number"
                    bind:value={r_peFargLenSd} />
                </div>
              {/if}
            {:else if r_readSimulator === 'badread'}
              <div class="field-row">
                <label class="field-label" for="r-brprofile" data-tooltip="Badread error model name. Examples: nanopore2023, nanopore2020, pacbio2016. Leave blank for the default model. See Badread documentation for the full list.">Sequencer profile / model <span class="tip">?</span></label>
                <input id="r-brprofile" class="text-input" type="text" bind:value={r_sequencerProfile}
                  placeholder="e.g. nanopore2023" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brcovdepth" data-tooltip="Target average sequencing depth (×) for Badread long-read simulation.">Coverage depth <span class="tip">?</span></label>
                <input id="r-brcovdepth" class="text-input short" type="number" step="0.1" min="0"
                  bind:value={r_coverageDepth} />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brlenm" data-tooltip="Mean length of simulated long reads (bp). Leave blank to use Badread's default distribution for the selected model.">Read length mean (bp) <span class="tip">?</span></label>
                <input id="r-brlenm" class="text-input short" type="number" min="1"
                  bind:value={r_longReadLenMean} placeholder="default" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brlens" data-tooltip="Standard deviation of long-read lengths (bp). Leave blank to use Badread defaults.">Read length SD (bp) <span class="tip">?</span></label>
                <input id="r-brlens" class="text-input short" type="number" min="0"
                  bind:value={r_longReadLenSd} placeholder="default" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brglitch" data-tooltip="Badread glitch parameters as 'rate,size,skips'. Simulates signal glitches in the sequencer. E.g. '10000,25,5' = 1 glitch per 10,000 bp, size 25, skipping 5 bases.">Glitches (rate,size,skips) <span class="tip">?</span></label>
                <input id="r-brglitch" class="text-input short" type="text"
                  bind:value={r_badreadGlitches} placeholder="e.g. 10000,25,5" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brjunk" data-tooltip="Percentage of reads that are entirely random sequence (junk). Simulates empty nanopores or other noise. Default: 0.">Junk reads (%) <span class="tip">?</span></label>
                <input id="r-brjunk" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_badreadJunkReads} placeholder="0" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brrandom" data-tooltip="Percentage of reads drawn from random sequence instead of the reference. Simulates background noise from sequencer. Default: 0.">Random reads (%) <span class="tip">?</span></label>
                <input id="r-brrandom" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_badreadRandomReads} placeholder="0" />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-brchim" data-tooltip="Percentage of chimeric reads (two random fragments joined end-to-end). Common in long-read datasets from adapter ligation artifacts. Default: 0.">Chimeric reads (%) <span class="tip">?</span></label>
                <input id="r-brchim" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_badreadChimeras} placeholder="0" />
              </div>
            {/if}
            <h4 class="subsection">Grouping (optional)</h4>
            <FilePicker label="Target groups" bind:value={r_groups}
              tooltip="Two-column TSV: sequence_id → group_name. Collapses multiple variant sequences of the same organism into one logical entity for TP/FP/FN classification."
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <FilePicker label="Distractor groups" bind:value={r_distractorGroups}
              tooltip="Two-column TSV: contig_id → group_name. By default, all contigs from each distractor FASTA are automatically grouped by filename. Use this file when multiple organisms are mixed into one FASTA."
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <h4 class="subsection">Host Filtering (optional)</h4>
            <FilePicker label="Host FASTA" bind:value={r_hostFasta}
              tooltip="Host genome FASTA (e.g. human). Reads aligning to this reference are removed before mapping to targets. Simulates host-depletion in a metagenomics workflow."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            {#if r_hostFasta}
              <div class="field-row">
                <label class="field-label" for="r-hostminimap" data-tooltip="Minimap2 preset for host filtering alignment. Auto-selected based on --minimap-preset by default. Override only if your host reads require a different preset.">Host minimap2 preset <span class="tip">?</span></label>
                <select id="r-hostminimap" class="text-input short" bind:value={r_hostMinimapPreset}>
                  <option value="sr">sr (Illumina short reads)</option>
                  <option value="map-ont">map-ont (Nanopore)</option>
                  <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                  <option value="asm5">asm5 (high-identity assembly)</option>
                  <option value="asm10">asm10</option>
                  <option value="asm20">asm20</option>
                  <option value="ava-ont">ava-ont (ONT overlap)</option>
                  <option value="ava-pb">ava-pb (PacBio overlap)</option>
                  <option value="cdna">cdna</option>
                  <option value="lr:hq">lr:hq</option>
                  <option value="lr:hqae">lr:hqae</option>
                  <option value="map-iclr">map-iclr</option>
                  <option value="map-pb">map-pb (PacBio CLR)</option>
                  <option value="splice">splice</option>
                  <option value="splice:hq">splice:hq</option>
                  <option value="splice:sr">splice:sr</option>
                </select>
              </div>
            {/if}
            <h4 class="subsection">Species Identification</h4>
            <label class="check-label">
              <input type="checkbox" bind:checked={r_identify} />
              Run species-level identification after metrics (genome mode)
            </label>
            {#if r_identify}
              <div class="field-row">
                <label class="field-label" for="r-identthresh" data-tooltip="Minimum sequence identity (%) between targets for cross-reactivity to be used when resolving ambiguous species calls.">Identity threshold (%) <span class="tip">?</span></label>
                <input id="r-identthresh" class="text-input short" type="number" min="0" max="100"
                  bind:value={r_identifyIdentThreshold} />
              </div>
              <div class="field-row">
                <label class="field-label" for="r-minuniq" data-tooltip="Minimum number of targets unique to a species that must be detected to confidently call it PRESENT rather than AMBIGUOUS.">Min unique targets for call <span class="tip">?</span></label>
                <input id="r-minuniq" class="text-input short" type="number" min="1"
                  bind:value={r_identifyMinUniq} />
              </div>
            {/if}
          </AdvancedOptions>

        <!-- ── build-probes ──────────────────────────── -->
        {:else if selectedTool === 'build-probes'}
          <section class="form-section" data-kind="input">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={bp_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section" data-kind="params">
            <h3>Design Method</h3>
            <div class="field-row">
              <label class="field-label" for="bp-method" data-tooltip="Tile: simple sliding window (fastest). CATCH-lite: greedy set cover targeting a coverage fraction. Syotti-lite: greedy set cover with k-mer seeding. ProbeTools-lite: iterative k-mer clustering. CATCH: external tool (requires catch conda package).">Method <span class="tip">?</span></label>
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
                <label class="field-label" for="bp-step" data-tooltip="Spacing between consecutive probes (bp). Positive = gap between probes. Negative = overlap. E.g. -20 = probes overlap by 20 bp.">Step (negative = overlap) <span class="tip">?</span></label>
                <input id="bp-step" class="text-input short" type="number" bind:value={bp_step} />
              </div>
            {:else if bp_method === 'catch-lite' || bp_method === 'catch'}
              <div class="field-row">
                <label class="field-label" for="bp-catchstride" data-tooltip="Step between candidate probe positions in the initial tiling before set-cover optimization. Smaller stride = more candidates but slower.">Probe stride <span class="tip">?</span></label>
                <input id="bp-catchstride" class="text-input short" type="number" bind:value={bp_catchStride} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchmm">Mismatches</label>
                <input id="bp-catchmm" class="text-input short" type="number" bind:value={bp_catchMismatches} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchext" data-tooltip="Extend each candidate probe by this many bp on each side before set-cover. Increases the genomic region each probe can cover.">Extension (bp each side) <span class="tip">?</span></label>
                <input id="bp-catchext" class="text-input short" type="number" min="0" bind:value={bp_catchExtension} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchcov" data-tooltip="Fraction of each target sequence that must be covered by at least one probe. Set-cover stops when this is achieved. 1.0 = full coverage required.">Min coverage fraction <span class="tip">?</span></label>
                <input id="bp-catchcov" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_catchCoverage} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchmh" data-tooltip="Probes with MinHash Jaccard similarity above this threshold are considered duplicates and merged during CATCH set-cover. Lower value = more aggressive deduplication.">MinHash dedup threshold <span class="tip">?</span></label>
                <input id="bp-catchmh" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_catchMinhashThreshold} />
              </div>
            {:else if bp_method === 'syotti-lite'}
              <div class="field-row">
                <label class="field-label" for="bp-syottimm" data-tooltip="Maximum mismatches allowed between a probe and a target region for the probe to count as covering it. Higher = more tolerant, fewer probes needed.">Mismatches <span class="tip">?</span></label>
                <input id="bp-syottimm" class="text-input short" type="number" bind:value={bp_syottiMismatches} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-syottiseed" data-tooltip="K-mer seed length for Syotti's greedy set-cover index. Shorter seeds find more hits but are slower; longer seeds are faster but may miss divergent sequences.">Seed length <span class="tip">?</span></label>
                <input id="bp-syottiseed" class="text-input short" type="number" bind:value={bp_syottiSeedLen} />
              </div>
            {:else if bp_method === 'probetools-lite'}
              <div class="field-row">
                <label class="field-label" for="bp-ptstep" data-tooltip="Step size for k-mer sampling when building the ProbeTools cluster index. Smaller step = denser sampling, better coverage, slower.">K-mer step <span class="tip">?</span></label>
                <input id="bp-ptstep" class="text-input short" type="number" min="1" bind:value={bp_ptStep} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptident" data-tooltip="Minimum sequence identity for two sequences to be assigned to the same ProbeTools cluster. 0.9 = 90% identity required to cluster together.">Cluster identity threshold <span class="tip">?</span></label>
                <input id="bp-ptident" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_ptIdentity} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptcov" data-tooltip="Fraction of each target sequence that must be covered by probes. Algorithm stops adding probes once this is reached per target.">Target coverage fraction <span class="tip">?</span></label>
                <input id="bp-ptcov" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_ptCoverage} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptbatch" data-tooltip="Number of probe candidates added per iteration. Larger batches are faster but may overshoot the optimal panel size.">Probes per iteration <span class="tip">?</span></label>
                <input id="bp-ptbatch" class="text-input short" type="number" min="1" bind:value={bp_ptBatchSize} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptmax" data-tooltip="Hard cap on the total number of probes in the final panel. Leave blank for no limit. Useful when you need to fit probes on a fixed-size array.">Max panel size (blank = unlimited) <span class="tip">?</span></label>
                <input id="bp-ptmax" class="text-input short" type="number" min="1" bind:value={bp_ptMaxPanelSize} placeholder="unlimited" />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptdepth" data-tooltip="Minimum number of probes that must cover each base for it to count as covered. Depth > 1 provides redundancy against probe failure.">Min coverage depth <span class="tip">?</span></label>
                <input id="bp-ptdepth" class="text-input short" type="number" min="1" bind:value={bp_ptMinDepth} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptiter" data-tooltip="Maximum number of probe-addition iterations before stopping, even if coverage goal is not met.">Max iterations <span class="tip">?</span></label>
                <input id="bp-ptiter" class="text-input short" type="number" min="1" bind:value={bp_ptMaxIterations} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-ptgain" data-tooltip="Stop iterating if the fractional coverage gain in the last iteration falls below this value. Prevents infinite loops when adding probes no longer helps.">Min coverage gain to continue <span class="tip">?</span></label>
                <input id="bp-ptgain" class="text-input short" type="number" step="0.0001" min="0" bind:value={bp_ptMinCoverageGain} />
              </div>
            {/if}
          </section>
          <section class="form-section" data-kind="output">
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
              <label class="field-label" for="bp-mingc" data-tooltip="Minimum GC content fraction for a probe to pass. Probes below this value are too AT-rich and may have low melting temperatures. Typical range: 0.30–0.70.">Min GC <span class="tip">?</span></label>
              <input id="bp-mingc" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_minGc} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-maxgc" data-tooltip="Maximum GC content fraction for a probe to pass. Probes above this value are too GC-rich and may form secondary structures or be hard to synthesize. Typical range: 0.30–0.70.">Max GC <span class="tip">?</span></label>
              <input id="bp-maxgc" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_maxGc} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-maxnfrac" data-tooltip="Target sequences with more than this fraction of ambiguous N bases are excluded before probe design. High-N sequences produce unreliable probes.">Max N fraction in targets <span class="tip">?</span></label>
              <input id="bp-maxnfrac" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_maxNFrac} />
            </div>
            <label class="check-label" data-tooltip="After design, replace each N in a probe with T (or A/C/G if adjacent to T). Avoids ambiguous bases that some synthesizers cannot manufacture.">
              <input type="checkbox" bind:checked={bp_noNInProbes} />
              Replace N bases in probes with non-N <span class="tip">?</span>
            </label>
            <div class="field-row">
              <label class="field-label" for="bp-dust" data-tooltip="sDUST complexity score threshold. Probes with a score above this are considered low-complexity (repetitive) and removed. Lower = stricter filtering. Default ~2.5.">sDUST threshold <span class="tip">?</span></label>
              <input id="bp-dust" class="text-input short" type="number" step="0.1" bind:value={bp_dustThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-dustwin" data-tooltip="Sliding window size for sDUST complexity scoring. Larger windows average over more sequence context; smaller windows flag short repetitive runs.">sDUST window size (bp) <span class="tip">?</span></label>
              <input id="bp-dustwin" class="text-input short" type="number" min="1" bind:value={bp_dustWindow} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-maxmask" data-tooltip="Maximum fraction of a probe that can be masked by sDUST before the probe is dropped. E.g. 0.5 = probes where more than half the sequence is low-complexity are removed.">Max masked fraction <span class="tip">?</span></label>
              <input id="bp-maxmask" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_maxMaskedFrac} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-collapse" data-tooltip="cd-hit-est identity threshold for collapsing redundant target sequences before probe design. Sequences above this identity are merged into a single representative. 0.95 = 95% identity.">Collapse identity threshold <span class="tip">?</span></label>
              <input id="bp-collapse" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_collapseThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-dedup" data-tooltip="cd-hit-est identity threshold for deduplicating final probes after design. Probes above this identity are collapsed to one representative. Removes near-identical probes that would waste array space.">Dedup identity threshold <span class="tip">?</span></label>
              <input id="bp-dedup" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={bp_dedupThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-thresh" data-tooltip="Minimum alignment identity (%) for a probe to be flagged as cross-reactive with a genome or another probe. Higher = only flag very close matches.">Xreact homology threshold (%) <span class="tip">?</span></label>
              <input id="bp-thresh" class="text-input short" type="number" min="0" max="100" bind:value={bp_threshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-aligner" data-tooltip="Aligner for the cross-reactivity step. minimap2 is fast and needs no external install, but its minimizer seeding can miss weak or short homologous regions. blast (blastn-short) is more sensitive, at the cost of speed — requires BLAST+ on PATH. Uses the Threads setting below.">Xreact aligner <span class="tip">?</span></label>
              <select id="bp-aligner" class="select-input" bind:value={bp_aligner}>
                <option value="minimap2">minimap2 (fast, embedded)</option>
                <option value="blast">blast (more sensitive, requires BLAST+)</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-prox" data-tooltip="Two cross-reactive hits within this distance (bp) on the same reference are merged into one report entry. Reduces noise from overlapping alignments.">Proximity distance (bp) <span class="tip">?</span></label>
              <input id="bp-prox" class="text-input short" type="number" min="0" bind:value={bp_proximity} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-minimap" data-tooltip="Minimap2 alignment preset for probe assessment. 'sr' works for most probe lengths (≤250 bp). Use 'asm5' for high-identity assemblies or 'map-ont' for long-read contexts.">Minimap2 preset <span class="tip">?</span></label>
              <select id="bp-minimap" class="text-input short" bind:value={bp_minimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
            </div>
            <h4 class="subsection">Coverage Refinement</h4>
            <div class="field-row">
              <label class="field-label" for="bp-refmode" data-tooltip="After initial probe design, optionally re-run assessment and prune probes that fall below the coverage threshold. 'iterations' = fixed rounds; 'stable' = repeat until no more probes are pruned.">Refinement <span class="tip">?</span></label>
              <select id="bp-refmode" class="select-input" bind:value={bp_refineMode}>
                <option value="none">None</option>
                <option value="iterations">Fixed iterations</option>
                <option value="stable">Until stable</option>
              </select>
            </div>
            {#if bp_refineMode !== 'none'}
              <div class="field-row">
                <label class="field-label" for="bp-refthresh" data-tooltip="A probe is pruned during refinement if removing it leaves at least this percentage of target bases still covered at 1×. Lower = more aggressive pruning.">1× coverage threshold (%) <span class="tip">?</span></label>
                <input id="bp-refthresh" class="text-input short" type="number" min="0" max="100" bind:value={bp_refineThreshold} />
              </div>
              {#if bp_refineMode === 'iterations'}
                <div class="field-row">
                  <label class="field-label" for="bp-refiter" data-tooltip="Number of refinement rounds to run. Each round re-assesses coverage and prunes redundant probes.">Iterations <span class="tip">?</span></label>
                  <input id="bp-refiter" class="text-input short" type="number" min="1" bind:value={bp_refineIterations} />
                </div>
              {/if}
            {/if}
            <h4 class="subsection">Assessment</h4>
            <FilePicker label="Genomes FASTA (cross-reactivity check)" bind:value={bp_genomes}
              tooltip="Optional: FASTA of non-target genomes (e.g. human, host) to check probes against for cross-reactivity. Probes hitting these are flagged in the assessment report."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <label class="check-label" data-tooltip="Skip the assess-probes step that runs after build-probes. Use if you only want the probe FASTA and will assess separately.">
              <input type="checkbox" bind:checked={bp_skipAssess} />
              Skip assessment step <span class="tip">?</span>
            </label>
            <div class="field-row">
              <label class="field-label" for="bp-threads">Threads</label>
              <input id="bp-threads" class="text-input short" type="number" min="1" bind:value={threads} />
            </div>
          </AdvancedOptions>

        <!-- ── assess-probes ────────────────────────── -->
        {:else if selectedTool === 'assess-probes'}
          <section class="form-section" data-kind="input">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={ap_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={ap_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Genomes FASTA (optional)" bind:value={ap_genomes}
              tooltip="Optional FASTA of off-target genomes (e.g. host) to include in cross-reactivity analysis. If omitted, only self cross-reactivity (probe-vs-probe) is checked."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section" data-kind="output">
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
              <label class="field-label" for="ap-threshold" data-tooltip="Minimum alignment identity (%) for a probe to be flagged as cross-reactive with a genome or another probe. Higher = only report very close matches.">Homology threshold (%) <span class="tip">?</span></label>
              <input id="ap-threshold" class="text-input short" type="number" min="0" max="100"
                bind:value={ap_threshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="ap-aligner" data-tooltip="Aligner for the cross-reactivity step. minimap2 is fast and needs no external install, but its minimizer seeding can miss weak or short homologous regions. blast (blastn-short) is more sensitive, at the cost of speed — requires BLAST+ on PATH. Uses the Threads setting below.">Xreact aligner <span class="tip">?</span></label>
              <select id="ap-aligner" class="select-input" bind:value={ap_aligner}>
                <option value="minimap2">minimap2 (fast, embedded)</option>
                <option value="blast">blast (more sensitive, requires BLAST+)</option>
              </select>
            </div>
            <label class="check-label" data-tooltip="By default, each target is assessed individually (probes aligned one target at a time) to distinguish true gaps from multi-mapper ambiguity. Skip this for panels with >10 000 targets to save time.">
              <input type="checkbox" bind:checked={ap_noIndividual} />
              Skip per-target individual coverage mapping (faster for large panels) <span class="tip">?</span>
            </label>
            <div class="field-row">
              <label class="field-label" for="ap-prox" data-tooltip="Cross-reactive hits within this distance (bp) on the same reference are merged into one entry. Reduces noise from overlapping alignments.">Proximity distance (bp) <span class="tip">?</span></label>
              <input id="ap-prox" class="text-input short" type="number" min="0" bind:value={ap_proximity} />
            </div>
            <div class="field-row">
              <label class="field-label" for="ap-minimap" data-tooltip="Minimap2 alignment preset. 'sr' works for probes ≤250 bp (default). Use 'map-ont' for long-read contexts or 'asm5' for near-identical assemblies.">Minimap2 preset <span class="tip">?</span></label>
              <select id="ap-minimap" class="text-input short" bind:value={ap_minimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="ap-gaplen" data-tooltip="Only report uncovered gaps in the detail TSV if they are at least this many base pairs long. Shorter gaps are too small to matter for most capture designs. Defaults to median probe length.">Min gap length in detail output (bp) <span class="tip">?</span></label>
              <input id="ap-gaplen" class="text-input short" type="number" min="1"
                bind:value={ap_gapMinLength} placeholder="default" />
            </div>
            <h4 class="subsection">Coverage Refinement</h4>
            <div class="field-row">
              <label class="field-label" for="ap-refmode" data-tooltip="Optionally re-assess and prune redundant probes: 'iterations' = fixed rounds; 'stable' = repeat until no more probes can be removed without dropping coverage below threshold.">Refinement <span class="tip">?</span></label>
              <select id="ap-refmode" class="select-input" bind:value={ap_refineMode}>
                <option value="none">None</option>
                <option value="iterations">Fixed iterations</option>
                <option value="stable">Until stable</option>
              </select>
            </div>
            {#if ap_refineMode !== 'none'}
              <div class="field-row">
                <label class="field-label" for="ap-refthresh" data-tooltip="A probe is pruned if removing it leaves at least this percentage of target bases still covered at 1×. Lower = more aggressive pruning of redundant probes.">1× coverage threshold (%) <span class="tip">?</span></label>
                <input id="ap-refthresh" class="text-input short" type="number" min="0" max="100" bind:value={ap_refineThreshold} />
              </div>
              {#if ap_refineMode === 'iterations'}
                <div class="field-row">
                  <label class="field-label" for="ap-refiter" data-tooltip="Number of refinement rounds. Each round removes probes that fall below the coverage threshold and reassesses.">Iterations <span class="tip">?</span></label>
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
          <section class="form-section" data-kind="input">
            <h3>Inputs</h3>
            <FilePicker label="Probes FASTA" bind:value={xr_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Against FASTA(s) (genomes to check against)" multiple
              bind:value={xr_against}
              tooltip="One or more FASTA files of non-target genomes to check probes against (e.g. human host, common contaminants). Probes with significant hits are flagged as cross-reactive."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <label class="check-label" data-tooltip="Also check probes against each other (probe-vs-probe). Flags probe pairs that may co-hybridize or compete for the same binding site.">
              <input type="checkbox" bind:checked={xr_self} />
              Probe-vs-probe self cross-reactivity <span class="tip">?</span>
            </label>
          </section>
          <section class="form-section" data-kind="output">
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
              <label class="field-label" for="xr-threshold" data-tooltip="Minimum alignment identity (%) to flag a probe as cross-reactive. 80% = flag alignments with ≥80% sequence identity. Higher values = more specific, fewer false flags.">Homology threshold (%) <span class="tip">?</span></label>
              <input id="xr-threshold" class="text-input short" type="number" min="0" max="100"
                bind:value={xr_threshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="xr-aligner" data-tooltip="minimap2 is fast and needs no external install, but its minimizer seeding can miss weak or short homologous regions. blast (blastn-short) is more sensitive to short/divergent homology, at the cost of speed — requires BLAST+ on PATH.">Aligner <span class="tip">?</span></label>
              <select id="xr-aligner" class="select-input" bind:value={xr_aligner}>
                <option value="minimap2">minimap2 (fast, embedded)</option>
                <option value="blast">blast (more sensitive, requires BLAST+)</option>
              </select>
            </div>
            {#if xr_aligner === 'blast'}
              <div class="field-row">
                <label class="field-label" for="xr-threads">Threads</label>
                <input id="xr-threads" class="text-input short" type="number" min="1" bind:value={threads} />
              </div>
            {:else}
            <div class="field-row">
              <label class="field-label" for="xr-minimap" data-tooltip="Minimap2 alignment preset for cross-reactivity alignment. 'sr' for short probes (≤250 bp), 'asm5' for near-identical sequences, 'map-ont' for long-read contexts.">Minimap2 preset <span class="tip">?</span></label>
              <select id="xr-minimap" class="text-input short" bind:value={xr_minimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
            </div>
            {/if}
          </AdvancedOptions>

        <!-- ── panel-qc ─────────────────────────────── -->
        {:else if selectedTool === 'panel-qc'}
          <section class="form-section" data-kind="input">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={pq_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map (TSV)" bind:value={pq_sampleTargetMap} required
              tooltip="TSV mapping species/genome IDs to their target sequence IDs (genome_id → target_id). Tells panel-qc which targets belong to the same species for discriminability scoring."
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
          </section>
          <section class="form-section" data-kind="output">
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
              <label class="field-label" for="pq-ident" data-tooltip="Minimum alignment identity (%) to count two targets as cross-similar. Used to compute inter-species confusion and discriminability scores. Higher = only flag very close pairs.">Identity threshold (%) <span class="tip">?</span></label>
              <input id="pq-ident" class="text-input short" type="number" min="0" max="100"
                bind:value={pq_identityThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="pq-minimap" data-tooltip="Minimap2 preset for target-vs-target similarity alignment. 'asm5' works well for comparing closely related viral or bacterial sequences.">Minimap2 preset <span class="tip">?</span></label>
              <select id="pq-minimap" class="text-input short" bind:value={pq_minimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
            </div>
          </AdvancedOptions>

        <!-- ── identify ─────────────────────────────── -->
        {:else if selectedTool === 'identify'}
          <section class="form-section" data-kind="input">
            <h3>Inputs</h3>
            <FilePicker label="detected_detail.tsv" bind:value={id_detectedDetail} required
              tooltip="The detected_detail.tsv output from a previous baitbench run (or baitbench metrics). Contains per-target read counts used for species calling."
              tooltipBelow
              filters={[{ name: 'TSV', extensions: ['tsv'] }]} />
            <FilePicker label="Sample-target map (TSV)" bind:value={id_sampleTargetMap} required
              tooltip="TSV mapping genome/species IDs to target sequence IDs. Tells identify which targets correspond to the same species for multi-target pattern calling."
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <FilePicker label="Pre-computed target similarity TSV (optional)"
              bind:value={id_targetSimilarity}
              tooltip="Pre-computed target-vs-target similarity matrix (from a prior panel-qc or assess-probes run). Speeds up identify by skipping the alignment step."
              filters={[{ name: 'TSV', extensions: ['tsv'] }]} />
            <FilePicker label="Targets FASTA (optional, compute similarity on-the-fly)"
              bind:value={id_targets}
              tooltip="Provide target sequences to compute cross-similarity on-the-fly (instead of a pre-computed TSV). Used to explain away false positives caused by cross-reactive probes."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section" data-kind="output">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="id-outprefix">Output prefix</label>
              <input id="id-outprefix" class="text-input" type="text" bind:value={outputPrefix} placeholder="optional" />
            </div>
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="id-ident" data-tooltip="Alignment identity threshold (%) for target similarity. Two species are considered cross-similar if their targets share this much identity — used to explain apparent false positives.">Identity threshold (%) <span class="tip">?</span></label>
              <input id="id-ident" class="text-input short" type="number" min="0" max="100"
                bind:value={id_identityThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="id-minuniq" data-tooltip="Minimum number of uniquely detected targets (not explained by cross-reactivity) required to call a species PRESENT. Higher = fewer false positives, but risks missing low-titer species.">Min unique targets for call <span class="tip">?</span></label>
              <input id="id-minuniq" class="text-input short" type="number" min="1"
                bind:value={id_minUniqueTargets} />
            </div>
            <div class="field-row">
              <label class="field-label" for="id-minimap" data-tooltip="Minimap2 preset for on-the-fly target similarity computation. 'asm5' is recommended for comparing closely related pathogen sequences.">Minimap2 preset <span class="tip">?</span></label>
              <select id="id-minimap" class="text-input short" bind:value={id_minimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
            </div>
          </AdvancedOptions>

        <!-- ── coverage-curve ────────────────────────── -->
        {:else if selectedTool === 'coverage-curve'}
          <section class="form-section" data-kind="input">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={cc_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={cc_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Distractor FASTA(s)" bind:value={cc_distractors} multiple required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Genomes FASTA (optional, genome mode)" bind:value={cc_genomes}
              tooltip="Full genome sequences (e.g. complete bacterial chromosomes). In genome mode, fragments are generated from these genomes but reads are mapped to targets. Leave blank for standard virus/amplicon mode."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map (optional)" bind:value={cc_sampleTargetMap}
              tooltip="TSV linking genome IDs to target sequence IDs (genome mode only). Enables correct read attribution when full genomes are simulated but reads map to gene targets."
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
          </section>

          <section class="form-section" data-kind="input">
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

          <section class="form-section" data-kind="params">
            <h3>Simulation</h3>
            <div class="field-row">
              <label class="field-label" for="cc-nfrags" data-tooltip="Total fragments to simulate per pipeline run. More fragments = more stable coverage estimates but slower runtime. 1 000–10 000 is typical.">Number of fragments <span class="tip">?</span></label>
              <input id="cc-nfrags" class="text-input short" type="number" min="100"
                bind:value={cc_numFragments} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-simmode" data-tooltip="Thermodynamic: uses SantaLucia nearest-neighbor model to bias fragment selection by probe binding free energy (ΔG). Simple: purely random fragment selection. Thermodynamic is more realistic but slower.">Simulate mode <span class="tip">?</span></label>
              <select id="cc-simmode" class="select-input" bind:value={cc_simulateMode}>
                <option value="thermodynamic">Thermodynamic (TNN)</option>
                <option value="simple">Simple</option>
              </select>
            </div>
          </section>

          <section class="form-section" data-kind="params">
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
                <label class="sweep-label" data-tooltip="qPCR cycle threshold values to sweep. Higher CT = lower target fraction in the mixture. Enter space-separated values (e.g. 20 25 30). Each value runs a separate pipeline.">
                  <input type="checkbox" bind:checked={cc_ctSweep} />
                  CT <span class="tip">?</span>
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
                    <label class="field-label" for="cc-ctb" data-tooltip="The reference CT value where target fraction equals the baseline fraction. Default: CT 20 = 1% target. Calibrate to your assay's known reference point.">CT baseline <span class="tip">?</span></label>
                    <input id="cc-ctb" class="text-input short" type="number" step="0.1" bind:value={cc_ctBaseline} />
                  </div>
                  <div class="field-row">
                    <label class="field-label" for="cc-ctbf" data-tooltip="Target fraction at the baseline CT. E.g. 0.01 means at CT 20, 1% of the mixture is target. The CT-to-fraction conversion is anchored to this point.">CT baseline fraction <span class="tip">?</span></label>
                    <input id="cc-ctbf" class="text-input short" type="number" step="0.001" min="0" max="1" bind:value={cc_ctBaselineFraction} />
                  </div>
                  <div class="field-row">
                    <label class="field-label" for="cc-cte" data-tooltip="PCR amplification efficiency per cycle. 1.0 = perfect doubling (100%). Real assays typically run at 0.90–0.98. Affects how steeply target fraction changes per CT unit.">PCR efficiency (0–1) <span class="tip">?</span></label>
                    <input id="cc-cte" class="text-input short" type="number" step="0.01" min="0" max="1" bind:value={cc_ctEfficiency} />
                  </div>
                {:else}
                  <p class="hint-sm">Provide two (CT, target-fraction) reference points; efficiency is derived automatically.</p>
                  <div class="field-row">
                    <label class="field-label" for="cc-ctcal1" data-tooltip="First calibration point as 'CT,fraction'. E.g. '20,0.01' means at CT 20, 1% of the mixture is target. PCR efficiency is derived from the two points automatically.">Point 1 (CT,fraction) <span class="tip">?</span></label>
                    <input id="cc-ctcal1" class="text-input short" type="text" bind:value={cc_ctCal1} placeholder="20,0.01" />
                  </div>
                  <div class="field-row">
                    <label class="field-label" for="cc-ctcal2" data-tooltip="Second calibration point as 'CT,fraction'. E.g. '25,0.003'. Together with Point 1, defines the CT-to-fraction conversion curve for your specific assay.">Point 2 (CT,fraction) <span class="tip">?</span></label>
                    <input id="cc-ctcal2" class="text-input short" type="text" bind:value={cc_ctCal2} placeholder="25,0.003" />
                  </div>
                {/if}
              </AdvancedOptions>
            {:else}
              <div class="sweep-row">
                <label class="sweep-label" data-tooltip="Fraction of fragments from distractor (off-target) sequences. 0.9 = 90% background noise. Enter space-separated values to sweep (e.g. 0.5 0.7 0.9).">
                  <input type="checkbox" bind:checked={cc_distractorFracSweep} />
                  Distractor Fraction <span class="tip">?</span>
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
                <label class="sweep-label" data-tooltip="Hybridization temperature for TNN thermodynamic scoring. Higher temp = more stringent probe selection. Enter space-separated values to sweep (e.g. 60 65 70 75).">
                  <input type="checkbox" bind:checked={cc_tempSweep} />
                  Hybridization Temp (°C) <span class="tip">?</span>
                </label>
                {#if cc_tempSweep}
                  <input class="text-input sweep-input" type="text" bind:value={cc_tempList} placeholder="60 65 70 75" />
                {:else}
                  <input class="text-input short" type="number" step="1" bind:value={cc_tempFixed} placeholder="70" />
                {/if}
              </div>
            {/if}

            <div class="sweep-row">
              <label class="sweep-label" data-tooltip="Fraction of probe-hybridizing fragments that are actually captured and sequenced. Models enrichment efficiency. 0.5 = 50% recovery. Sweep to compare sensitivity across efficiency levels.">
                <input type="checkbox" bind:checked={cc_cfSweep} />
                Capture Fraction <span class="tip">?</span>
              </label>
              {#if cc_cfSweep}
                <input class="text-input sweep-input" type="text" bind:value={cc_cfList} placeholder="0.3 0.5 0.8" />
              {:else}
                <input class="text-input short" type="number" min="0" max="1" step="0.01" bind:value={cc_cfFixed} placeholder="0.5" />
              {/if}
            </div>

            <div class="sweep-row">
              <label class="sweep-label" data-tooltip="Number of reads/sequences to subsample before mapping. Simulates different sequencing depths. Leave blank (or uncheck) to use all simulated reads. Sweep to plot depth-vs-coverage curves.">
                <input type="checkbox" bind:checked={cc_nsSweep} />
                Num Sequences <span class="tip">?</span>
              </label>
              {#if cc_nsSweep}
                <input class="text-input sweep-input" type="text" bind:value={cc_nsList} placeholder="500 1000 5000" />
              {:else}
                <input class="text-input short" type="number" min="1" bind:value={cc_nsFixed} placeholder="all" />
              {/if}
            </div>
          </section>

          <section class="form-section" data-kind="output">
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
              <label class="field-label" for="cc-seed" data-tooltip="Fixed seed for the random number generator. Set to the same value to reproduce identical results across runs. Leave blank for a different result each time.">Random seed (blank = random) <span class="tip">?</span></label>
              <input id="cc-seed" class="text-input short" type="text" bind:value={cc_seed} placeholder="e.g. 42" />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-fraglenmean" data-tooltip="Mean length of simulated library fragments (before sequencing). Should match your expected library insert size. Typical WGS: 300–500 bp; capture: 200–400 bp.">Fragment length mean (bp) <span class="tip">?</span></label>
              <input id="cc-fraglenmean" class="text-input short" type="number" bind:value={cc_fragLenMean} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-fraglenmin" data-tooltip="Minimum allowed fragment length. Fragments shorter than this are rejected during sampling. Should be ≥ probe length for meaningful capture simulation.">Fragment length min (bp) <span class="tip">?</span></label>
              <input id="cc-fraglenmin" class="text-input short" type="number" bind:value={cc_fragLenMin} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-fraglenmax" data-tooltip="Maximum allowed fragment length. Very long fragments are less efficiently captured; capping here keeps the simulation realistic for your size-selection range.">Fragment length max (bp) <span class="tip">?</span></label>
              <input id="cc-fraglenmax" class="text-input short" type="number" bind:value={cc_fragLenMax} />
            </div>
            <FilePicker label="Host FASTA (optional filtering)" bind:value={cc_hostFasta}
              tooltip="Optional host genome for pre-mapping host read removal. Reads that align here are discarded before mapping to targets. Use to simulate depletion of host background."
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <div class="field-row">
              <label class="field-label" for="cc-minimap" data-tooltip="Minimap2 preset for mapping reads to targets. 'sr' for short reads (default), 'map-ont' for long reads, 'asm5' for high-identity assemblies.">Minimap2 preset <span class="tip">?</span></label>
              <select id="cc-minimap" class="text-input short" bind:value={cc_minimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-hostminimap" data-tooltip="Separate minimap2 preset for host filtering alignment. Override only if host reads require a different preset than the main mapping.">Host minimap2 preset <span class="tip">?</span></label>
              <select id="cc-hostminimap" class="text-input short" bind:value={cc_hostMinimapPreset}>
                <option value="sr">sr (Illumina short reads)</option>
                <option value="map-ont">map-ont (Nanopore)</option>
                <option value="map-hifi">map-hifi (PacBio HiFi)</option>
                <option value="asm5">asm5 (high-identity assembly)</option>
                <option value="asm10">asm10</option>
                <option value="asm20">asm20</option>
                <option value="ava-ont">ava-ont (ONT overlap)</option>
                <option value="ava-pb">ava-pb (PacBio overlap)</option>
                <option value="cdna">cdna</option>
                <option value="lr:hq">lr:hq</option>
                <option value="lr:hqae">lr:hqae</option>
                <option value="map-iclr">map-iclr</option>
                <option value="map-pb">map-pb (PacBio CLR)</option>
                <option value="splice">splice</option>
                <option value="splice:hq">splice:hq</option>
                <option value="splice:sr">splice:sr</option>
              </select>
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
    padding: 0 18px;
    height: 46px;
    background: #1e3a8a;
    flex-shrink: 0;
    position: relative;
    z-index: 10;
  }
  .header-side {
    flex: 1;
    display: flex;
    align-items: center;
  }
  .header-side.right { justify-content: flex-end; }
  .logo {
    font-size: 1.4rem;
    font-weight: 800;
    color: #ffffff;
    letter-spacing: -0.02em;
    user-select: none;
  }
  .btn-running {
    padding: 4px 10px;
    background: rgba(255, 255, 255, 0.15);
    border: 1px solid rgba(255, 255, 255, 0.3);
    border-radius: 5px;
    font-size: 0.8rem;
    color: #ffffff;
    cursor: pointer;
    white-space: nowrap;
    transition: background 0.15s;
  }
  .btn-running:hover { background: rgba(255, 255, 255, 0.22); }
  .body {
    display: flex;
    flex: 1;
    overflow: hidden;
  }
  /* Sidebar */
  .sidebar {
    width: 188px;
    flex-shrink: 0;
    border-right: none;
    background: var(--sidebar-bg);
    padding: 6px 0;
    overflow-y: auto;
    display: flex;
    flex-direction: column;
    gap: 0;
  }
  .cat-label {
    font-size: 0.65rem;
    font-weight: 700;
    text-transform: uppercase;
    color: var(--sidebar-cat-text);
    letter-spacing: 0.09em;
    padding: 16px 14px 5px;
    border-top: 1px solid var(--sidebar-separator);
    margin-top: 6px;
  }
  .sidebar > :global(div:first-child) {
    border-top: none;
    margin-top: 0;
    padding-top: 8px;
  }
  .tool-btn {
    width: 100%;
    text-align: left;
    padding: 7px 18px;
    font-size: 0.85rem;
    border: none;
    background: transparent;
    cursor: pointer;
    color: var(--sidebar-text);
    border-radius: 0;
    transition: background 0.1s;
  }
  .tool-btn:hover { background: var(--sidebar-hover); }
  .tool-btn.active {
    background: var(--sidebar-active-bg);
    color: var(--sidebar-active-text);
    font-weight: 600;
  }
  .sidebar-spacer { flex: 1; }
  .sidebar-utilities {
    padding: 8px 12px 14px;
    border-top: 1px solid var(--sidebar-separator);
    display: flex;
    flex-direction: column;
    gap: 5px;
  }
  .util-btn {
    width: 100%;
    text-align: left;
    padding: 7px 12px;
    font-size: 0.8rem;
    border: 1px solid rgba(255, 255, 255, 0.18);
    background: rgba(255, 255, 255, 0.07);
    border-radius: 6px;
    cursor: pointer;
    color: var(--sidebar-text);
    transition: background 0.15s, border-color 0.15s;
    box-sizing: border-box;
  }
  .util-btn:hover {
    background: rgba(255, 255, 255, 0.13);
    border-color: rgba(255, 255, 255, 0.3);
  }
  /* Form area */
  .form-area {
    flex: 1;
    display: flex;
    flex-direction: column;
    overflow: hidden;
  }
  .tool-header {
    padding: 16px 22px 12px;
    border-bottom: 2px solid var(--color-border);
    background: var(--color-card);
    flex-shrink: 0;
  }
  .tool-header h2 { margin: 0 0 4px; font-size: 1.15rem; font-weight: 700; }
  .tool-desc { margin: 0; font-size: 0.82rem; color: var(--color-muted); }
  .form-scroll {
    flex: 1;
    overflow-y: auto;
    padding: 14px 16px 28px;
    display: flex;
    flex-direction: column;
    gap: 10px;
  }
  .form-section {
    display: flex;
    flex-direction: column;
    gap: 10px;
    background: var(--color-card);
    border: 1px solid var(--color-border);
    border-radius: 8px;
    padding: 14px 16px;
  }
  .form-section h3 {
    margin: 0 0 4px;
    font-size: 0.72rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.07em;
    color: var(--color-label);
    padding-bottom: 8px;
    border-bottom: 1px solid var(--color-border);
  }
  .form-section[data-kind="input"] {
    border-left: 3px solid rgba(59, 130, 246, 0.65);
    background: rgba(59, 130, 246, 0.03);
  }
  .form-section[data-kind="params"] {
    border-left: 3px solid rgba(139, 92, 246, 0.55);
    background: rgba(139, 92, 246, 0.03);
  }
  .form-section[data-kind="output"] {
    border-left: 3px solid rgba(16, 185, 129, 0.65);
    background: rgba(16, 185, 129, 0.03);
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
</style>
