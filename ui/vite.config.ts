import { defineConfig } from "vite";

// The workbench is a static single-page app. In development it serves from the
// repository root's `public` equivalent; the workbench.json artifact is written
// by `diachron export-workbench` and placed alongside the built assets.
export default defineConfig({
  root: ".",
  build: {
    outDir: "dist",
    emptyOutDir: true,
  },
  server: {
    port: 5173,
  },
});
