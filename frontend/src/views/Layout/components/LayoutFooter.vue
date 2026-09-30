<template>
  <footer class="site-footer">
    <div class="container footer-inner">
      <div class="footer-brand">
        <strong>SCAID</strong>
        <span>Single-Cell AutoImmune Disease database · read-only, no account required</span>
      </div>
      <nav class="footer-links" aria-label="Footer">
        <router-link to="/Information">Information</router-link>
        <a href="/api/schema/">API</a>
        <router-link to="/TeamInformation">Atlases</router-link>
        <a href="https://github.com/c-chunhao/SCAID-portal" target="_blank" rel="noopener noreferrer">Source code</a>
        <a href="mailto:chunhao.bio@gmail.com">chunhao.bio@gmail.com</a>
      </nav>
    </div>
  </footer>
</template>

<style scoped lang="scss">
.site-footer {
  background: var(--scaid-soft, #f8f5f2);
  border-top: 1px solid var(--scaid-line, #e5ddd7);
  color: var(--scaid-muted, #6b5f5a);
  font-size: .8125rem;
  line-height: 1.6;
  width: 100%;
}

.footer-inner {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px 28px;
  flex-wrap: wrap;
  padding-top: 24px;
  padding-bottom: 24px;
}

.footer-brand {
  display: flex;
  gap: 6px 12px;
  align-items: baseline;
  flex-wrap: wrap;
  min-width: 0;
}

.footer-brand strong {
  color: var(--scaid-accent, #8a503a);
  letter-spacing: .04em;
}

.footer-brand span { max-width: 62ch; }

.footer-links {
  display: flex;
  align-items: center;
  gap: 0 20px;
  flex-wrap: wrap;
}

.footer-links a {
  display: inline-flex;
  align-items: center;
  min-height: 44px;
  color: var(--scaid-ink, #2b2320);
  text-decoration: none;
  text-underline-offset: 4px;
  overflow-wrap: anywhere;
}

.footer-links a:hover {
  color: var(--scaid-accent, #8a503a);
  text-decoration: underline;
}

.footer-links a:focus-visible {
  outline: 2px solid var(--scaid-accent, #8a503a);
  outline-offset: 3px;
  border-radius: 2px;
}

@media (max-width: 575px) {
  .footer-inner { gap: 12px; }
  .footer-brand { flex-direction: column; gap: 4px; }
}
</style>
