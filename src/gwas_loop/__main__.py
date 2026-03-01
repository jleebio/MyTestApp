"""CLI entry point: python -m gwas_loop"""
import argparse
import logging
from .runner import AutonomousRunner, LoopConfig


def main():
    parser = argparse.ArgumentParser(description="GWAS Self-Improving Research Loop")
    parser.add_argument("--traits", nargs="+", default=None, help="Trait IDs to process")
    parser.add_argument("--max-iter", type=int, default=10, help="Max loop iterations")
    parser.add_argument("--data-dir", default="data", help="Data directory")
    parser.add_argument("--output-dir", default="reports", help="Output directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate pipeline only")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    config = LoopConfig(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        max_iterations=args.max_iter,
    )

    runner = AutonomousRunner(config)
    print(f"Pipeline status: {runner.status()}")
    results = runner.run(traits=args.traits, dry_run=args.dry_run)
    print(f"\nCompleted {len(results)} iterations")
    for r in results:
        print(f"  [{r.trait}] iter={r.iteration}: {r.n_hypotheses_generated} hypotheses")


if __name__ == "__main__":
    main()
