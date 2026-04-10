import { describe, expect, it, vi } from "vitest";
import {
  compileRegex,
  expandWithNeighbors,
  findMatchingNodeIds,
  getNodeText,
  waitFor,
  type VisNetwork,
} from "./graphSearch";

describe("graphSearch", () => {
  it("getNodeText concatenates id/label/title", () => {
    expect(
      getNodeText({ id: "X", label: "L", title: "<b>T</b>" }).trim().split("\n")
    ).toEqual(["X", "L", "<b>T</b>"]);
  });

  it("compileRegex throws on empty", () => {
    expect(() => compileRegex("   ")).toThrow(/empty/i);
  });

  it("findMatchingNodeIds matches id/label/title", () => {
    const nodes = [
      { id: "A", label: "alpha", title: "hello" },
      { id: "B", label: "beta", title: "world" },
      { id: "C", label: "gamma", title: "Oligo: weak" },
    ];
    // Search matches anywhere in the combined text, like the UI does.
    const re = compileRegex("beta");
    expect(findMatchingNodeIds(nodes, re)).toEqual(["B"]);
    expect(findMatchingNodeIds(nodes, compileRegex("oligo:\\s*weak"))).toEqual(["C"]);
  });

  it("expandWithNeighbors includes 1-hop neighbors and de-dupes", () => {
    const network: VisNetwork = {
      getConnectedNodes: (id) => (id === "A" ? ["B", "C"] : id === "B" ? ["A"] : []),
      selectNodes: () => undefined,
      focus: () => undefined,
      fit: () => undefined,
      unselectAll: () => undefined,
    };
    const ids = expandWithNeighbors(network, ["A", "B"]);
    expect(new Set(ids)).toEqual(new Set(["A", "B", "C"]));
  });

  it("waitFor resolves when value becomes available", async () => {
    const seq = [null, null, "ok"] as Array<string | null>;
    const resolve = vi.fn(() => seq.shift());
    const v = await waitFor(resolve, { timeoutMs: 1000, intervalMs: 1 });
    expect(v).toBe("ok");
    expect(resolve).toHaveBeenCalled();
  });

  it("waitFor times out", async () => {
    await expect(waitFor(() => null, { timeoutMs: 5, intervalMs: 1 })).rejects.toThrow(/timeout/i);
  });
});

