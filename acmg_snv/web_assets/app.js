"use strict";

const byId = (id) => document.getElementById(id);
const form = byId("analysis");
const classes = {
  P: "Classe 5 · Pathogène",
  LP: "Classe 4 · Probablement pathogène",
  VUS: "Classe 3 · Signification incertaine",
  LB: "Classe 2 · Probablement bénin",
  B: "Classe 1 · Bénin",
};
const strengths = {
  supporting: "En soutien",
  moderate: "Modérée",
  strong: "Forte",
  very_strong: "Très forte",
  standalone: "Indépendante (BA1)",
};
const criterionCodes = [
  "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5",
  "PM6", "PP1", "PP2", "PP3", "PP4", "BA1", "BS1", "BS2", "BS3", "BS4",
  "BP1", "BP2", "BP4", "BP5", "BP7",
];
const reportStatuses = {
  needs_context: "Contexte à compléter",
  out_of_scope: "Hors périmètre",
  not_evaluable: "Classe non évaluable",
  provisional: "Classification provisoire",
  review_required: "Relecture requise",
  approved: "Classification relue",
};
const criterionStatuses = {
  applied: "Appliqué",
  removed: "Retiré",
  excluded: "Exclu",
  overridden: "Modifié manuellement",
  not_applied: "Non appliqué",
  not_met: "Non rempli",
  not_evaluable: "Non évaluable",
  conflict: "En conflit",
  manual: "Appliqué manuellement",
  suppressed: "Écarté",
  not_applicable: "Non applicable",
};
const relationships = {
  unknown: "Phase inconnue",
  trans: "En trans",
  cis: "En cis",
  homozygous: "Homozygote",
};
const phaseChecks = {
  affected: "Cas atteint avec phénotype concordant",
  rare: "SNV sous le seuil PM2 de l’affection",
  partner_rare: "Autre allèle sous le seuil PM2 de l’affection",
  independent_classification: "Classe du partenaire établie indépendamment",
  phase_confirmed: "Phase ou homozygotie confirmée",
  same_condition: "Partenaire dans le même gène pour la même affection",
};
const contextFields = [
  "inheritance", "gene", "disease", "mechanism", "phenotype", "segregation",
  "transcript", "conflict-resolution", "vcep-exclusion-reason",
];
const contextLabels = {
  gene: "Gène", transcript: "Transcrit", disease: "Affection",
  inheritance: "Transmission", mechanism: "Mécanisme",
  phenotype: "Phénotype", phase: "Phase",
  segregation: "Ségrégation", target_genes: "Gènes étudiés",
  conflict_resolution: "Résolution des conflits",
  vcep_exclusion_reason: "Justification de l’exclusion VCEP",
};
let functionalTsv = "";
let functionalReady = true;
let currentReport = null;
let pending = false;
let editVersion = 0;

function element(tag, text, className) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (className) node.className = className;
  return node;
}

function select(options, label) {
  const control = element("select");
  if (label) control.add(new Option(label, ""));
  for (const [value, text] of Object.entries(options)) {
    control.add(new Option(text, value));
  }
  return control;
}

function addField(container, title, key, control) {
  const label = element("label", title);
  control.dataset.key = key;
  label.append(control);
  container.append(label);
  return control;
}

function addCheck(container, title, key) {
  const label = element("label", undefined, "checkbox");
  const control = element("input");
  control.type = "checkbox";
  control.dataset.key = key;
  label.append(control, document.createTextNode(title));
  container.append(label);
}

function addText(container, title, key, placeholder = "") {
  const control = element("input");
  control.placeholder = placeholder;
  return addField(container, title, key, control);
}

function addRemove(container, parent = container) {
  const button = element("button", "Retirer", "secondary");
  button.type = "button";
  button.addEventListener("click", () => {
    parent.remove();
    invalidateReview();
  });
  container.append(button);
}

function manualFields(row, item) {
  const grid = element("div", undefined, "row");
  const codes = Object.fromEntries(criterionCodes.map((code) => [code, code]));
  const actions = {apply: "Appliquer ou modifier", remove: "Écarter"};
  const code = addField(grid, "Critère", "code", select(codes, "Choisir"));
  const strength = select(strengths, "Sélectionner");
  addField(grid, "Force", "strength", strength);
  const action = addField(grid, "Décision", "action", select(actions));
  addRemove(grid, row);
  row.append(grid);
  if (item.code) code.value = item.code;
  if (item.strength in strengths) strength.value = item.strength;
  action.addEventListener("change", () => {
    strength.disabled = action.value === "remove";
  });
}

function manualRow(item = {}) {
  const row = element("fieldset", undefined, "manual-row");
  manualFields(row, item);
  addText(row, "Justification", "reason");
  addText(row, "Source", "source");
  addText(row, "Identifiants des preuves", "evidence_ids", "Virgules");
  byId("manual-rows").append(row);
  invalidateReview();
  return row;
}

function addPhaseChecks(row) {
  const details = element("details");
  details.append(element("summary", "Éléments pour évaluer la phase"));
  const checks = element("div", undefined, "field-grid");
  for (const [key, label] of Object.entries(phaseChecks)) {
    addCheck(checks, label, key);
  }
  details.append(checks);
  const options = select({yes: "Oui", no: "Non"}, "Non vérifiée");
  addField(details, "Cooccurrence surtout en cis", "predominantly_cis",
    options);
  addText(details, "Identifiant local du cas", "case_id", "Facultatif");
  row.append(details);
}

function phaseRow() {
  const row = element("fieldset", undefined, "phase-row");
  const grid = element("div", undefined, "row");
  const labels = Object.fromEntries(
    Object.entries(classes).map(([code, label]) => [code, label.slice(11)]),
  );
  addText(grid, "Variant associé", "variant", "Coordonnées ou HGVS");
  addField(grid, "Relation", "relationship", select(relationships));
  const classification = select(labels, "Non précisée");
  addField(grid, "Classe du partenaire", "classification", classification);
  addRemove(grid, row);
  row.append(grid);
  addText(row, "Justification", "reason");
  addPhaseChecks(row);
  byId("phase-rows").append(row);
  invalidateReview();
  row.querySelector("input").focus();
}

function fieldValues(row, selector = "[data-key]") {
  const entries = Array.from(row.querySelectorAll(selector), (field) => {
    let value = field.value.trim();
    if (field.type === "checkbox") value = field.checked;
    if (field.type === "number" && value !== "") value = field.valueAsNumber;
    return [field.dataset.key, value];
  });
  return Object.fromEntries(entries.filter(([, value]) => value !== ""));
}

function rowValues(selector) {
  return Array.from(document.querySelectorAll(selector), (row) => {
    return fieldValues(row);
  });
}

function clinicalPanel(title, id) {
  const panel = element("details", undefined, "clinical-panel");
  panel.id = id;
  panel.append(element("summary", title));
  byId("clinical-structured").append(panel);
  return panel;
}

function addEntry(panel, label, build) {
  const list = element("div");
  const button = element("button", label, "secondary");
  button.type = "button";
  button.addEventListener("click", () => {
    const row = build();
    addRemove(row);
    list.append(row);
    invalidateReview();
    row.querySelector("input").focus();
  });
  panel.append(list, button);
}

function buildPhenotypePanel() {
  const panel = clinicalPanel("Preuve phénotypique · PP4", "pp4-evidence");
  const specificity = select({
    specific: "Spécifique de l’affection",
    additional_features: "Signes supplémentaires discriminants",
    pathognomonic: "Pathognomonique",
  }, "À préciser");
  addField(panel, "Spécificité du phénotype", "specificity", specificity);
  addCheck(panel, "Principales causes de ce phénotype explorées",
    "majority_causes_tested");
  addCheck(panel, "Spécificité élevée confirmée en RCP", "mdt_confirmed");
  addText(panel, "Justification PP4", "reason");
  addText(panel, "Identifiant local du cas", "case_id", "Facultatif");
}

function segregationFamily() {
  const row = element("fieldset", undefined, "segregation-family");
  addText(row, "Identifiant local de la famille", "family_id");
  addText(row, "Probabilité de coségrégation", "probability", "1/8");
  const meioses = element("input");
  meioses.type = "number";
  meioses.min = 1;
  addField(row, "Ou nombre de méioses informatives (AD)",
    "informative_meioses", meioses);
  addCheck(row, "Hypothèses du calcul simplifié en AD vérifiées",
    "jarvik_simple_assumptions");
  return row;
}

function buildSegregationPanel() {
  const panel = clinicalPanel("Coségrégation · PP1", "pp1-evidence");
  addCheck(panel, "Lien gène-maladie définitif", "definitive_gene_disease");
  addCheck(panel, "Coségrégation documentée", "cosegregation");
  addCheck(panel, "Familles indépendantes", "independent_families");
  addText(panel, "Justification PP1", "reason");
  addEntry(panel, "Ajouter une famille", segregationFamily);
}

function denovoCase() {
  const row = element("fieldset", undefined, "denovo-case");
  addText(row, "Identifiant local du cas", "case_id", "Facultatif");
  const parents = {confirmed: "Confirmés", assumed: "Présumés"};
  const parentage = select(parents, "À préciser");
  addField(row, "Liens de parenté", "parental_relationships", parentage);
  const fit = select({specific: "Spécifique", consistent: "Concordant",
    heterogeneous: "Concordant et génétiquement hétérogène",
    inconsistent: "Discordant"}, "À préciser");
  addField(row, "Concordance du phénotype", "phenotype_fit", fit);
  addCheck(row, "Les deux parents ont été testés", "parents_tested");
  addCheck(row, "Apparition de novo établie", "denovo_established");
  addCheck(row, "Variant P ou LP supplémentaire en AR", "additional_plp");
  addCheck(row, "Mosaïcisme germinal", "germline_mosaicism");
  addText(row, "Justification de la preuve de novo", "reason");
  addText(row, "Source", "source");
  return row;
}

function buildDenovoPanel() {
  const panel = clinicalPanel("Preuve de novo · PS2 / PM6", "denovo-evidence");
  addCheck(panel, "Cas indépendants", "denovo_unrelated");
  addEntry(panel, "Ajouter un cas", denovoCase);
}

function clinicalValues() {
  const phenotype = fieldValues(byId("pp4-evidence"));
  const directControls = ":scope > label > input";
  const segregation = fieldValues(byId("pp1-evidence"), directControls);
  const independent = byId("denovo-evidence").querySelector("input").checked;
  segregation.families = rowValues(".segregation-family");
  const denovo = rowValues(".denovo-case");
  const clinical = {};
  if (denovo.length) clinical.denovo = denovo;
  if (independent) clinical.denovo_unrelated = true;
  if (phenotype.reason || phenotype.specificity) {
    clinical.phenotype_evidence = phenotype;
  }
  if (segregation.reason || segregation.families.length) {
    clinical.segregation_evidence = segregation;
  }
  return clinical;
}

function splitEntries(text) {
  return text.split(/[,\n;]+/).map((value) => value.trim()).filter(Boolean);
}

function manualValues() {
  return rowValues(".manual-row").map((item) => {
    const apply = item.action === "apply";
    const missingApplication = apply && (!item.strength || !item.source);
    if (!item.code || !item.reason || missingApplication) {
      byId("manual-section").open = true;
      const message = "Justifiez la décision et renseignez force et source.";
      throw new Error(message);
    }
    return {...item, evidence_ids: splitEntries(item.evidence_ids || "")};
  });
}

function phaseValues() {
  const values = rowValues(".phase-row");
  const incomplete = ({variant, relationship}) => {
    return !variant && relationship !== "homozygous";
  };
  if (values.some(incomplete)) {
    byId("clinical-section").open = true;
    throw new Error("Identifiez chaque variant associé ou retirez sa ligne.");
  }
  return values.map(({predominantly_cis: mostlyCis, ...item}) => {
    if (mostlyCis) item.predominantly_cis = mostlyCis === "yes";
    return item;
  });
}

function advancedContext() {
  const text = byId("advanced-context").value.trim();
  if (!text) return {};
  try {
    const context = JSON.parse(text);
    if (context && typeof context === "object" && !Array.isArray(context)) {
      return context;
    }
  } catch {
    // Le même message couvre la syntaxe JSON et le type de la valeur.
  }
  byId("clinical-section").open = true;
  throw new Error("Le contexte avancé doit être un objet JSON valide.");
}

function contextValues() {
  const entries = contextFields.map((id) => [
    id.replaceAll("-", "_"), byId(id).value.trim(),
  ]);
  const simple = Object.fromEntries(entries.filter(([, value]) => value));
  const context = {...advancedContext(), ...simple, ...clinicalValues()};
  const targetGenes = splitEntries(byId("targeted-genes").value);
  const phases = phaseValues();
  if (targetGenes.length) context.target_genes = targetGenes;
  if (phases.length) context.phase = phases;
  if (byId("phase-unrelated").checked) context.phase_unrelated = true;
  context.reviewed = byId("reviewed").checked;
  context.vcep_reviewed = byId("vcep-reviewed").checked;
  context.vcep_population_confirmed = byId("vcep-population-confirmed").checked;
  const specification = byId("vcep-specification").value;
  if (specification) context.vcep_specification = specification;
  return context;
}

function getPayload() {
  return {
    variant: byId("variant").value.trim(),
    context: contextValues(),
    annotation: byId("annotation").value.trim(),
    overrides: manualValues(),
    ps3: functionalTsv,
    ps3_name: byId("functional-file").files[0]?.name || "",
  };
}

async function loadFunctional(event) {
  const file = event.target.files[0];
  const status = byId("functional-status");
  functionalTsv = "";
  functionalReady = !file;
  status.textContent = "";
  byId("clear-functional").hidden = !file;
  if (!file) return;
  if (file.size > 24 * 1024 * 1024) {
    status.textContent = "Choisissez un tableau de 24 Mio maximum.";
    return;
  }
  try {
    const bytes = await file.arrayBuffer();
    if (event.target.files[0] !== file) return;
    const decoder = new TextDecoder("utf-8", {fatal: true});
    functionalTsv = decoder.decode(bytes);
    functionalReady = true;
    status.textContent = `${file.name} prêt pour l’analyse locale.`;
  } catch {
    if (event.target.files[0] === file) {
      status.textContent = "Ce fichier n’a pas pu être lu en UTF-8.";
    }
  }
}

function clearFunctional() {
  const input = byId("functional-file");
  input.value = "";
  loadFunctional({target: input});
  invalidateReview();
}

function updateContextPrompt(report) {
  const prompt = byId("context-prompt");
  const inheritanceMissing = !report.context?.inheritance;
  prompt.hidden = report.status !== "needs_context" || !inheritanceMissing;
  if (!prompt.hidden) {
    prompt.textContent = "Sélectionnez le mode de transmission pour évaluer.";
    byId("inheritance").focus();
  }
  if (report.context?.gene && !byId("gene").value) {
    byId("gene").placeholder = `Gène proposé : ${report.context.gene}`;
  }
  if (report.context?.transcript && !byId("transcript").value) {
    byId("transcript").placeholder = `Proposé : ${report.context.transcript}`;
  }
  const matches = report.vcep?.matches || [];
  const exactMatch = matches.some((match) => match.exact);
  byId("disease-prompt").hidden = !matches.length || exactMatch;
}

function appendList(container, title, values) {
  if (!values || !values.length) return;
  const list = element("ul");
  for (const item of values) list.append(element("li", readable(item)));
  container.append(element("h3", title), list);
}

function readable(value) {
  if (value === null || value === undefined || value === "") return "—";
  if (Array.isArray(value)) return value.map(readable).join("\n") || "—";
  return typeof value === "object" ? JSON.stringify(value) : String(value);
}

function reportTable(headers) {
  const wrapper = element("div", undefined, "table-scroll");
  const table = element("table");
  const row = table.createTHead().insertRow();
  headers.forEach((heading) => row.append(element("th", heading)));
  wrapper.append(table);
  return {wrapper, body: table.createTBody()};
}

function addValues(table, values) {
  const row = table.body.insertRow();
  for (const value of values) {
    row.append(element("td", readable(value), "plain-value"));
  }
  return row;
}

function editFromReport(item) {
  byId("manual-section").open = true;
  const existing = Array.from(document.querySelectorAll(".manual-row")).find(
    (row) => row.querySelector('[data-key="code"]').value === item.code,
  );
  const row = existing || manualRow(item);
  row.querySelector('[data-key="reason"]').focus();
}

function renderCriteria(container, criteria) {
  if (!criteria || !criteria.length) return;
  const table = reportTable([
    "Critère", "Évaluation", "Force", "Points", "Justification et sources",
    "",
  ]);
  for (const item of criteria) {
    const row = addValues(table, [
      item.code, criterionStatuses[item.status] || item.status,
      strengths[item.strength] || item.strength,
      item.points ?? (item.standalone ? "Hors score" : null),
      [item.reason, item.source, ...(item.evidence_ids || [])].filter(Boolean),
    ]);
    const cell = row.insertCell();
    if (criterionCodes.includes(item.code)) {
      const edit = element("button", "Modifier", "secondary report-edit");
      edit.type = "button";
      edit.setAttribute("aria-label", `Modifier le critère ${item.code}`);
      edit.addEventListener("click", () => editFromReport(item));
      cell.append(edit);
    }
  }
  container.append(element("h3", "Évaluation des critères"), table.wrapper);
}

function renderVcepCriteria(panel, criteria) {
  if (!criteria) return;
  const details = element("details");
  details.append(element("summary", "Spécifications par critère"));
  if (typeof criteria === "string") {
    details.append(element("p", criteria, "plain-value"));
  } else {
    details.append(element("pre", JSON.stringify(criteria, null, 2)));
  }
  panel.append(details);
}

function renderVcep(container, vcep) {
  updateVcepOptions(vcep?.matches || []);
  const guidance = vcep?.guidance || [];
  if (!guidance.length) return;
  const panel = element("details", undefined, "card");
  panel.append(element("summary", "Spécifications ClinGen VCEP"));
  for (const match of guidance) {
    const name = [match.panel, match.version, match.id].filter(Boolean);
    panel.append(element("h3", name.join(" · ") || "Spécification"));
    if (match.source) {
      panel.append(element("p", readable(match.source), "plain-value"));
    }
    renderVcepCriteria(panel, match.criteria);
  }
  container.append(panel);
}

function updateVcepOptions(matches) {
  const control = byId("vcep-specification");
  const previous = control.value;
  const keys = new Set();
  control.length = 1;
  for (const match of matches) {
    const key = `${match.id}@${match.version}`;
    if (keys.has(key)) continue;
    keys.add(key);
    control.add(new Option(`${key} · ${match.panel}`, key));
  }
  control.value = keys.has(previous) ? previous : "";
}

function renderCandidates(container, report) {
  if (report.status !== "needs_context" || !report.candidates?.length) return;
  const table = reportTable(["Gène", "Transcrit", "MANE", "Conséquence", ""]);
  for (const item of report.candidates) {
    const row = addValues(table, [
      item.gene, item.transcript, item.mane, item.consequence,
    ]);
    const choose = element("button", "Choisir", "secondary report-edit");
    choose.type = "button";
    choose.addEventListener("click", () => {
      byId("gene").value = item.gene || "";
      byId("transcript").value = item.transcript || "";
      invalidateReview();
      byId("evaluate").focus();
    });
    row.insertCell().append(choose);
  }
  container.append(element("h3", "Annotations candidates"), table.wrapper);
}

function renderContext(container, report) {
  const table = reportTable(["Élément", "Valeur"]);
  const variant = report.variant;
  if (variant) {
    const {chrom, pos, ref, alt} = variant;
    const display = `chr${chrom}:${pos}:${ref}:${alt}`;
    const coords = variant.key || display;
    addValues(table, ["Variant GRCh38", coords]);
  }
  for (const [key, label] of Object.entries(contextLabels)) {
    const value = report.context?.[key];
    if (value?.length || (value && typeof value !== "string")) {
      addValues(table, [label, value]);
    }
  }
  const reference = report.reference_checked ? "Vérifiée" : "—";
  addValues(table, ["Référence FASTA", reference]);
  addValues(table, ["Référentiel", report.framework]);
  addValues(table, ["Version du logiciel", report.software_version]);
  container.append(table.wrapper);
}

function renderReport(report) {
  currentReport = report;
  const container = byId("result-body");
  const category = report.classification;
  const classification = classes[category] ?? category;
  const status = reportStatuses[report.status] || report.status;
  const score = Number.isFinite(report.score) ? `Score : ${report.score}` : "";
  const summary = [classification, score, status].filter(Boolean).join(" · ");
  container.replaceChildren(element("p", summary, "result-summary"));
  if (report.score_reason) container.append(element("p", report.score_reason));
  renderContext(container, report);
  appendList(container, "Points à vérifier", report.alerts);
  appendList(container, "Conflits", report.conflicts);
  renderCandidates(container, report);
  renderVcep(container, report.vcep);
  renderCriteria(container, report.criteria);
  renderAnnotation(container, report.annotation);
  byId("output").hidden = false;
  byId("outdated").hidden = true;
  byId("download").disabled = false;
  byId("print").disabled = false;
}

function renderAnnotation(container, annotation) {
  if (!annotation) return;
  const details = element("details", undefined, "annotation-report");
  details.append(element("summary", "Annotation retenue"));
  details.append(element("pre", JSON.stringify(annotation, null, 2)));
  container.append(details);
}

function markOutdated() {
  editVersion += 1;
  if (!currentReport) return;
  byId("outdated").hidden = false;
  byId("download").disabled = true;
  byId("print").disabled = true;
}

function invalidateReview() {
  byId("reviewed").checked = false;
  byId("vcep-reviewed").checked = false;
  byId("vcep-population-confirmed").checked = false;
  markOutdated();
}

function onEdit(event) {
  if (event.target.closest("#review-section")) markOutdated();
  else invalidateReview();
}

async function requestClassification() {
  const response = await fetch("/classify", {
    method: "POST",
    headers: {"Content-Type": "application/json"},
    body: JSON.stringify(getPayload()),
    credentials: "omit",
  });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || "La requête a échoué.");
  return data;
}

async function evaluate(event) {
  event.preventDefault();
  if (pending) return;
  if (!functionalReady) {
    showStatus("Attendez la lecture du PS3 ou retirez le fichier.", true);
    return;
  }
  const version = editVersion;
  pending = true;
  byId("evaluate").disabled = true;
  showStatus("Évaluation en cours…");
  try {
    const report = await requestClassification();
    if (version === editVersion) showCompletedReport(report);
    else showStatus("Saisie modifiée. Relancez l’évaluation.");
  } catch (error) {
    showStatus(error.message, true);
  } finally {
    pending = false;
    byId("evaluate").disabled = false;
  }
}

function showCompletedReport(report) {
  renderReport(report);
  updateContextPrompt(report);
  if (byId("context-prompt").hidden) byId("result-title").focus();
  showStatus("");
}

function showStatus(message, error = false) {
  const status = byId("status");
  status.className = error ? "error" : "";
  status.textContent = message;
}

function downloadReport() {
  if (!currentReport || !byId("outdated").hidden) return;
  const contents = JSON.stringify(currentReport, null, 2) + "\n";
  const blob = new Blob([contents], {type: "application/json"});
  const url = URL.createObjectURL(blob);
  const download = element("a");
  download.href = url;
  download.download = "rapport_acmg.json";
  download.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

buildPhenotypePanel();
buildSegregationPanel();
buildDenovoPanel();
form.addEventListener("submit", evaluate);
form.addEventListener("input", onEdit);
byId("inheritance").addEventListener("change", () => {
  byId("context-prompt").hidden = true;
});
byId("functional-file").addEventListener("change", loadFunctional);
byId("clear-functional").addEventListener("click", clearFunctional);
byId("add-phase").addEventListener("click", phaseRow);
byId("add-criterion").addEventListener("click", () => {
  manualRow().querySelector("select").focus();
});
byId("download").addEventListener("click", downloadReport);
byId("print").addEventListener("click", () => window.print());
