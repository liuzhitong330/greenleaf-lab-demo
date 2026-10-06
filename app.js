(function () {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const numeric = (n) => typeof n === 'number' && Number.isFinite(n);
  const decimal = (n, places = 2) => numeric(n) ? n.toFixed(places) : 'not estimable';
  const percent = (n) => numeric(n) ? `${(100 * n).toFixed(1)}%` : 'not estimable';
  const pp = (n) => numeric(n) ? `${(100 * n).toFixed(1)} pp` : 'not estimable';
  const escape = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const number = (n) => numeric(n) ? n.toLocaleString('en-US') : '—';
  const svgNS = 'http://www.w3.org/2000/svg';
  let data, guide, target, showAll = false;

  function reasons(t) {
    const floor = Number($('read-floor').value);
    const tolerance = Number($('control-limit').value);
    const r = [];
    if (!numeric(t.baseline_reads) || t.baseline_reads < floor) r.push('low baseline support');
    if (!numeric(t.control_shift)) r.push('control comparison not estimable');
    else if (t.control_shift > tolerance) r.push('control-sensitive');
    if (!numeric(t.baseline_shift)) r.push('baseline comparison not estimable');
    else if (t.baseline_shift > tolerance) r.push('baseline-sensitive');
    if (!numeric(t.late_pair_gap)) r.push('late change not estimable');
    else if (t.late_pair_gap > tolerance) r.push('late change: inspect timing / extend coverage');
    return r;
  }

  function element(name, attrs, text) {
    const e = document.createElementNS(svgNS, name);
    Object.entries(attrs || {}).forEach(([k, v]) => e.setAttribute(k, v));
    if (text !== undefined) e.textContent = text;
    return e;
  }

  function bounds(values, fallback = [0, 1]) {
    const v = values.filter(numeric);
    let low = Math.min(fallback[0], ...v), high = Math.max(fallback[1], ...v);
    if (low === high) high = low + 1;
    const pad = (high - low) * 0.06;
    return [low - pad, high + pad];
  }

  function frame(svg, title, description, xb, yb, xlabel, ylabel, height = 330) {
    svg.replaceChildren(element('title', { id: `${svg.id}-title` }, title), element('desc', { id: `${svg.id}-desc` }, description));
    svg.setAttribute('aria-labelledby', `${svg.id}-title ${svg.id}-desc`);
    const plot = { left: 65, right: 582, top: 32, bottom: height - 53 };
    const x = (v) => plot.left + (v - xb[0]) / (xb[1] - xb[0]) * (plot.right - plot.left);
    const y = (v) => plot.bottom - (v - yb[0]) / (yb[1] - yb[0]) * (plot.bottom - plot.top);
    for (let i = 0; i <= 4; i++) {
      const v = yb[0] + i / 4 * (yb[1] - yb[0]);
      svg.append(element('line', { x1: plot.left, x2: plot.right, y1: y(v), y2: y(v), stroke: '#e5e5e5' }));
      svg.append(element('text', { x: plot.left - 9, y: y(v) + 4, 'text-anchor': 'end', 'font-size': 11 }, decimal(v, Math.abs(v) >= 100 ? 0 : 2)));
    }
    for (let i = 0; i <= 4; i++) {
      const v = xb[0] + i / 4 * (xb[1] - xb[0]);
      svg.append(element('text', { x: x(v), y: plot.bottom + 19, 'text-anchor': 'middle', 'font-size': 11 }, decimal(v, Number.isInteger(v) ? 0 : 1)));
    }
    svg.append(element('line', { x1: plot.left, x2: plot.right, y1: plot.bottom, y2: plot.bottom, stroke: '#999' }));
    if (yb[0] < 0 && yb[1] > 0) svg.append(element('line', { x1: plot.left, x2: plot.right, y1: y(0), y2: y(0), stroke: '#aaa', 'stroke-dasharray': '2 4' }));
    svg.append(element('text', { x: (plot.left + plot.right) / 2, y: height - 9, 'text-anchor': 'middle', 'font-size': 12 }, xlabel));
    svg.append(element('text', { x: plot.left, y: 14, 'font-size': 11 }, ylabel));
    return { x, y, plot };
  }

  function drawTimecourse() {
    const svg = $('hero-viz');
    const sets = [
      { values: target.f_cm, color: '#777', dash: '5 4', name: 'CM only' },
      { values: target.f_ko, color: '#c07a2a', dash: '1 5', name: 'KO only' },
      { values: target.f_bound, color: '#1f7a8c', dash: '', name: 'CM + KO' }
    ];
    const { x, y } = frame(svg, `${guide.name} ${target.mutation}: association measurements`, 'Unclipped apparent fraction depleted under pooled, CM-only and KO-only control normalization. Points are sequenced samples; they are not independent biological replicates.', [0, Math.max(...guide.times)], bounds(sets.flatMap((s) => s.values)), 'Actual collection time (minutes)', 'Fraction depleted / inferred binding (dimensionless)');
    sets.forEach((s) => {
      let path = '', started = false;
      s.values.forEach((v, i) => {
        if (!numeric(v)) { started = false; return; }
        path += `${started ? 'L' : 'M'}${x(guide.times[i]).toFixed(2)},${y(v).toFixed(2)} `;
        started = true;
      });
      svg.append(element('path', { d: path, fill: 'none', stroke: s.color, 'stroke-width': s.name === 'CM + KO' ? 2 : 1.4, 'stroke-dasharray': s.dash }));
      s.values.forEach((v, i) => {
        if (!numeric(v)) return;
        const p = element('circle', { cx: x(guide.times[i]), cy: y(v), r: s.name === 'CM + KO' ? 3 : 2, fill: s.color });
        p.append(element('title', {}, `${s.name}; ${guide.times[i]} min; sample ${guide.timepoint_ids[i]}; depleted ${decimal(v, 4)}; target ${target.counts[i]} reads; CM ${guide.controls.CM[i]}; KO ${guide.controls.KO[i]}`));
        svg.append(p);
      });
    });
    const flagged = reasons(target);
    $('readout').innerHTML = `<strong>${escape(guide.name)} · ${escape(target.mutation)}:</strong> late inferred binding ${percent(target.late)}; CM–KO difference ${pp(target.control_shift)}; baseline leave-one-out spread ${pp(target.baseline_shift)}; late change ${pp(target.late_pair_gap)}.<br>${flagged.length ? `Review: ${escape(flagged.join('; '))}.` : 'No review criterion exceeded at the current settings. This does not establish assay validity.'}`;
    $('sequence-note').innerHTML = `Library ${escape(guide.id)} · target ${escape(target.id)} · sequence includes the PAM<br><code>${escape(target.sequence)}</code>`;
    $('metric-one').textContent = number(target.baseline_reads);
    $('metric-two').textContent = percent(target.late);
    $('data-caption').innerHTML = `${guide.times.length} recorded samples; ${guide.times.filter((t) => t === 0).length} zero-time samples shown at their actual shared time. Late = mean of the last two samples (${guide.times.slice(-2).join(' and ')} min). Source: <a href="https://doi.org/10.6084/m9.figshare.12526157.v2">Boyle et al., Figshare v2</a>. Counts are reads, not biological replicates. Lines connect measurements, not fitted kinetics.`;
    renderRawTable();
  }

  function renderRawTable() {
    const tbody = $('raw-body');
    if (!tbody) return;
    tbody.replaceChildren();
    guide.times.forEach((time, i) => {
      const tr = document.createElement('tr');
      [guide.timepoint_ids[i], time, target.counts[i], guide.controls.CM[i], guide.controls.KO[i], decimal(target.f_bound[i], 4)].forEach((v) => {
        const td = document.createElement('td'); td.textContent = v; tr.append(td);
      });
      tbody.append(tr);
    });
  }

  function chooseTarget(id, scroll) {
    const selected = guide.targets.find((t) => t.id === id);
    if (!selected) return;
    target = selected;
    $('target-select').value = target.id;
    render();
    if (scroll) $('target-select').focus({ preventScroll: true });
  }

  function renderQueue() {
    const rows = guide.targets.map((t) => ({ target: t, reasons: reasons(t) }));
    const flagged = rows.filter((r) => r.reasons.length);
    flagged.sort((a, b) => b.reasons.length - a.reasons.length || (b.target.control_shift || 0) - (a.target.control_shift || 0));
    const displayed = showAll ? rows : flagged.slice(0, 8);
    const global = data.guides.reduce((total, g) => total + g.targets.filter((t) => reasons(t).length).length, 0);
    $('metric-three').textContent = `${flagged.length} / ${guide.targets.length}`;
    $('queue-summary').textContent = `${flagged.length} of ${guide.targets.length} targets enter this library’s review queue; ${number(global)} of ${number(data.meta.coverage.drilldown_targets)} across all ${data.guides.length} sublibraries at these settings. ${showAll ? 'All targets shown.' : `Showing ${Math.min(8, flagged.length)} priority entries.`} Review is not rejection: continued late change can reflect real kinetics.`;
    const tbody = $('queue-body'); tbody.replaceChildren();
    displayed.forEach((r) => {
      const tr = document.createElement('tr');
      if (r.target.id === target.id) tr.className = 'selected-row';
      const first = document.createElement('td');
      const button = document.createElement('button'); button.type = 'button'; button.className = 'text-button'; button.textContent = r.target.mutation;
      button.setAttribute('aria-label', `Inspect ${r.target.mutation}, target ${r.target.id}`);
      button.addEventListener('click', () => chooseTarget(r.target.id, true)); first.append(button); tr.append(first);
      [number(r.target.baseline_reads), pp(r.target.control_shift), r.reasons.join('; ') || 'No criterion exceeded'].forEach((v) => { const td = document.createElement('td'); td.textContent = v; tr.append(td); });
      tbody.append(tr);
    });
    if (!displayed.length) {
      const tr = document.createElement('tr'), td = document.createElement('td'); td.colSpan = 4; td.textContent = 'No targets in this library cross the current review criteria. View all targets or compare another library.'; tr.append(td); tbody.append(tr);
    }
    $('show-all').textContent = showAll ? 'Show review shortlist' : `Show all ${guide.targets.length} targets`;
    $('show-all').setAttribute('aria-expanded', String(showAll));
  }

  function drawModel() {
    const candidates = guide.targets.filter((t) => t.mutation !== 'WT' && numeric(t.model_ddg) && numeric(t.late));
    const svg = $('model-viz');
    const { x, y } = frame(svg, `${guide.name}: model penalty versus measured depletion`, 'Each point is one single-mismatch target. Selecting a point loads its raw time course. This same-study model comparison is not held-out validation.', bounds(candidates.map((t) => t.model_ddg), [0, 1]), bounds(candidates.map((t) => t.late)), 'Casland productive-binding penalty (kT)', 'Late fraction depleted / inferred binding', 320);
    candidates.forEach((t) => {
      const selected = t.id === target.id;
      const p = element('circle', { cx: x(t.model_ddg), cy: y(t.late), r: selected ? 5 : 3.5, fill: selected ? '#1f7a8c' : '#aaa', 'fill-opacity': selected ? 1 : 0.65, stroke: reasons(t).length ? '#444' : 'none', 'stroke-width': 1, tabindex: 0, role: 'button', class: 'chart-point', 'aria-label': `${t.mutation}: penalty ${decimal(t.model_ddg)} kT, late depleted ${percent(t.late)}. Inspect time course.` });
      p.append(element('title', {}, `${t.mutation}; ${decimal(t.model_ddg, 3)} kT; ${percent(t.late)} depleted${reasons(t).length ? '; review criterion exceeded' : ''}`));
      p.addEventListener('click', () => chooseTarget(t.id, true));
      p.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); chooseTarget(t.id, true); } });
      svg.append(p);
    });
    $('model-caption').innerHTML = `${candidates.length} single-mismatch targets from ${escape(guide.id)}. A dark outline means a current review criterion is exceeded, not a failed measurement. Select a point or use the target dropdown above. Relative energy penalties are mapped from author Casland parameters, with the mapping checked against <a href="https://github.com/augustboyle/Casland/tree/9193bf679c2bab346d8c0f30a1f92e52c29c9bed">the pinned code</a>; they are not absolute binding probabilities.`;
    renderChallenge(candidates);
  }

  function renderChallenge(candidates) {
    const wellSupported = candidates.filter((t) => !reasons(t).length).sort((a, b) => a.model_ddg - b.model_ddg || a.id.localeCompare(b.id));
    const review = candidates.filter((t) => reasons(t).length).sort((a, b) => (b.control_shift || 0) - (a.control_shift || 0));
    const items = [];
    const wt = guide.targets.find((t) => t.mutation === 'WT');
    if (wt) items.push([wt, 'Perfect-target reference', 'Anchor the sequence comparison; retain its own count and control checks.']);
    if (wellSupported.length) {
      const low = wellSupported[0], high = wellSupported[wellSupported.length - 1];
      if (high.model_ddg > low.model_ddg) {
        items.push([low, 'Lower model penalty', 'One edge of the model-penalty range among targets without a current review flag.']);
        items.push([high, 'Higher model penalty', 'The other edge of that range; a prospective challenge, not a predicted acceptance result.']);
      } else items.push([low, 'Supported model-penalty reference', 'The currently supported targets share one model penalty; they cannot establish a low–high challenge contrast.']);
    }
    if (review.length) items.push([review[0], 'Diagnostic review case', `${reasons(review[0]).join('; ')}. Recheck before giving the contrast a biological interpretation.`]);
    const box = $('challenge-panel'); box.replaceChildren();
    const intro = document.createElement('p'); intro.textContent = `Proposed ${items.length}-target panel for ${guide.name}. Selection updates with your review thresholds. Add the CM and KO normalization controls separately; they are not included in this target count.`; box.append(intro);
    const list = document.createElement('ol');
    items.forEach(([t, title, why]) => {
      const li = document.createElement('li'), button = document.createElement('button'); button.type = 'button'; button.className = 'text-button'; button.textContent = `${title}: ${t.mutation}`;
      button.addEventListener('click', () => chooseTarget(t.id, true)); li.append(button, document.createTextNode(` — ${decimal(t.model_ddg)} kT; ${percent(t.late)} late depletion. ${why}`)); list.append(li);
    });
    box.append(list);
    if (!wellSupported.length) { const p = document.createElement('p'); p.textContent = 'No single-mismatch target clears these descriptive review settings. Resolve the flagged measurements before using this library to span a prospective challenge panel.'; box.append(p); }
    else if (wellSupported[0].model_ddg === wellSupported[wellSupported.length - 1].model_ddg) { const p = document.createElement('p'); p.textContent = 'The supported targets do not span a model-penalty range at these settings. Resolve the relevant review signals or add independently supported sequences before testing a low–high contrast; do not relax thresholds solely to manufacture one.'; box.append(p); }
  }

  function render() { drawTimecourse(); renderQueue(); drawModel(); }

  function chooseGuide(id) {
    guide = data.guides.find((g) => g.id === id) || data.guides[0];
    $('guide-select').value = guide.id;
    $('target-select').replaceChildren();
    guide.targets.forEach((t) => { const o = document.createElement('option'); o.value = t.id; o.textContent = `${t.mutation === 'WT' ? 'Perfect target (WT)' : t.mutation} · ${t.id}`; $('target-select').append(o); });
    target = guide.targets.find((t) => t.mutation === 'WT') || guide.targets[0]; showAll = false; render();
  }

  function downloadQueue() {
    const fields = ['library', 'guide', 'target', 'mutation', 'sequence', 'baseline_reads', 'late_fraction_depleted', 'control_difference_fraction', 'baseline_leave_one_out_range', 'late_two_sample_change', 'casland_penalty_kT', 'baseline_floor', 'tolerance_fraction', 'review_reasons'];
    const quote = (v) => `"${String(v === null || v === undefined ? '' : v).replace(/"/g, '""')}"`;
    const rows = guide.targets.map((t) => [guide.id, guide.name, t.id, t.mutation, t.sequence, t.baseline_reads, t.late, t.control_shift, t.baseline_shift, t.late_pair_gap, t.model_ddg, $('read-floor').value, $('control-limit').value, reasons(t).join('; ') || 'No criterion exceeded']);
    const blob = new Blob([fields.map(quote).join(',') + '\r\n' + rows.map((r) => r.map(quote).join(',')).join('\r\n')], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob), a = document.createElement('a'); a.href = url; a.download = `greenleaf-${guide.id}-review-${$('read-floor').value}-${$('control-limit').value}.csv`; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function init(d) {
    if (!d || !Array.isArray(d.guides) || !d.guides.length) throw new Error('The published data package is unavailable.');
    data = d;
    data.guides.forEach((g) => { const o = document.createElement('option'); o.value = g.id; o.textContent = `${g.name} · ${g.id}`; $('guide-select').append(o); });
    const suggestions = [['S60', 'L1 reference'], ['S36', 'L1_CT comparison'], ['S18', 'Control-sensitive (S18)'], ['S1', 'VEGFA'], ['S32', 'Two baseline samples']];
    const label = document.createElement('span'); label.textContent = 'Compare:'; $('suggestions').append(label);
    suggestions.forEach(([id, name]) => { if (!data.guides.some((g) => g.id === id)) return; const b = document.createElement('button'); b.type = 'button'; b.textContent = name; b.addEventListener('click', () => chooseGuide(id)); $('suggestions').append(b); });
    $('guide-select').addEventListener('change', () => chooseGuide($('guide-select').value));
    $('target-select').addEventListener('change', () => chooseTarget($('target-select').value));
    ['read-floor', 'control-limit'].forEach((id) => $(id).addEventListener('change', render));
    $('show-all').addEventListener('click', () => { showAll = !showAll; renderQueue(); });
    $('download-queue').addEventListener('click', downloadQueue);
    $('normalization-note').textContent = `Coverage: ${number(data.meta.coverage.drilldown_targets)} target–library curves and ${number(data.meta.coverage.raw_count_measurements)} raw count measurements across ${data.guides.length} sublibraries. For each library, the baseline target/control ratio pools its metadata-listed zero-time samples (three for 90 libraries; two for S32). Expected target reads at time t = pooled baseline target/control ratio × control reads at t. Fraction depleted = 1 − observed target reads / (expected target reads + 0.1). CM denotes the complemented control and KO the PAM-disrupted control. The primary estimate pools CM + KO; alternatives use CM or KO alone. Control sensitivity is the absolute difference between the two alternative late estimates. Baseline sensitivity is the range of late estimates after leaving each baseline sample out. Late change is the absolute difference between the last two estimates; it may indicate ongoing kinetics. Estimates are never clipped.`;
    chooseGuide('S36');
  }

  Promise.resolve(window.GREENLEAF_DATA_READY || window.GREENLEAF_DATA).then(init).catch((error) => {
    $('readout').textContent = `The data could not be loaded. Please reload the public page or inspect the source repository. ${error.message}`;
    ['guide-select', 'target-select', 'read-floor', 'control-limit', 'show-all', 'download-queue'].forEach((id) => { $(id).disabled = true; });
  });
})();
