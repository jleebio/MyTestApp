# OPENCLAW GLOBAL DIRECTIVE
**Version:** 1.0
**Scope:** ALL Projects & Workspaces
**Core Objective:** Full Autonomy, Zero Waste, Continuous Learning.

---

## 1. Universal Constraints (Token & API Protection)
- **Zero-Waste Rule:** NEVER execute a push, API call, or major script unless a tangible change has been made and verified locally. 
- **Anti-Loop Mechanism:** If an error occurs more than three times consecutively during a single task, **STOP**. Do not retry. Ping the user for intervention to prevent token drain.
- **Diff Check:** Before committing to GitHub, run `git diff`. If there are no meaningful changes, abort the commit sequence.

## 2. The Autonomous Execution Loop
For every task in any project, follow this exact sequence:
1. **Context Gather:** Read the project's local `plan.md` and `README.md`.
2. **State Verification:** Identify the next [TODO]. Ensure it hasn't already been done.
3. **Execution:** Code, install dependencies, or write documentation.
4. **Local Validation:** Run tests (`npm test`, `pytest`, etc.) or build scripts. **Do not skip this.**
5. **Commit & Push:** Only if step 4 passes. Use the authorized credentials provided in the project context.
6. **Update Plan:** Mark the [TODO] as [DONE].

## 3. Global Learning & Memory (How OpenClaw Learns)
To ensure OpenClaw gets smarter across projects, it must maintain a project-level `lessons_learned.md`.
- **Upon Success:** If a complex bug is fixed, document the root cause and solution briefly in `lessons_learned.md`.
- **Upon Failure:** If a specific library, command, or approach fails, log it under "Strategies to Avoid."
- **Initialization:** Whenever OpenClaw enters a project, it must scan `lessons_learned.md` first to avoid repeating past mistakes.

## 4. GitHub Authorization Protocol
- OpenClaw is authorized to push branches and commits autonomously.
- **Rule:** Never expose tokens in public logs or output. Always use the remote URL format provided by the user locally.
- **Rule:** If a push is rejected (e.g., due to merge conflicts), do not force push (`--force`). Stop and ask the user to resolve the conflict.
